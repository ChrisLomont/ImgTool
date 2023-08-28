#pragma once
#include <array>
#include <memory>
#include <numbers> // c++20, for pi
#include <functional>
#include <cmath>
#include <vector>

#include "Image.h"

/* todo
 * reuse prev filter defs, integrate them all?
 * template the resampler on filter type, make it internally?
 */

namespace Lomont::Resampler2 {

    struct FilterDef
    {
        [[nodiscard]] virtual double support() const { return 0.0; }
        [[nodiscard]] virtual double filter(double t) const = 0;
        virtual ~FilterDef() = default;
    };
    struct BoxFilter : FilterDef
    {
        [[nodiscard]] double support() const override { return 0.5; }
        [[nodiscard]] double filter(double t) const override
        {
            if (-0.5 < t && t <= 0.5) return 1.0;
            return 0.0;
        }
    };
    struct BilinearFilter : FilterDef
    {
        [[nodiscard]] double support() const override { return 1.0; }
        [[nodiscard]] double filter(double t) const override
        {
            t = abs(t);
            if (t < 1.0) return 1.0 - t;
            return 0.0;
        }
    };
    struct BicubicFilter : FilterDef
    { // todo - B,C parametrize
        std::array<double, 4> c1{};
        std::array<double, 4> c2{};
        // Mitchell–Netravali parametrization (also called BC-splines)
        // B,C values:
        // 0, any -> Cardinal splines
        // 0, 1/2 -> Catmull-Rom spline, used in GIMP
        // 0, 3/4 -> used in photoshop
        // 1/3, 1/3 -> Mitchell-Netravali, used in ImageMagick
        // 1, 0 -> B-Spline, used in Paint.NET
        BicubicFilter(double B, double C)
        {
            c1 = { 6 - 2 * B, 0,-18 + 12 * B + 6 * C,  12 - 9 * B - 6 * C };
            c2 = { 8 * B + 24 * C, -12 * B - 48 * C, 6 * B + 30 * C, -B - 6 * C };
            for (auto i = 0; i < 4; ++i)
            { // normalize
                c1[i] /= 6;
                c2[i] /= 6;
            }
        }
        // Keys parametrization:
        // this is setting B=0 in MN parametrization, C = a
        BicubicFilter(double a = -0.5)
        {
            c1 = { 1,0,-a - 3,a + 2 };
            c2 = { -4 * a,8 * a,-5 * a,a };
        }

        [[nodiscard]] double support() const override { return 2.0; }
        [[nodiscard]] double filter(double t) const override
        {
            t = abs(t);
            const double t2 = t * t, t3 = t2 * t;
            if (t < 1.0)
                return c1[3] * t3 + c1[2] * t2 + c1[1] * t + c1[0];
            if (t < 2.0)
                return c2[3] * t3 + c2[2] * t2 + c2[1] * t + c2[0];
            return 0.0;
        }
    };
    struct LanczosFilter : FilterDef
    {
        double a{ 3 };
        LanczosFilter(double a = 3) : a(a) { }
        [[nodiscard]] double support() const override { return a; }
        [[nodiscard]] double filter(double t) const override
        {
            t = std::abs(t);
            if (a <= t) return 0.0;
            if (t < 1e-20) return 1.0;
            const auto pt = t * std::numbers::pi;
            if (t < 1e-10)
            {
                const auto a2 = a * a;
                // taylor series
                // 1 - (((1 + a^2) Pi^2) t^2)/(6 a ^ 2) + ((3 + 10 a ^ 2 + 3 a ^ 4) Pi ^ 4 t ^ 4) / (360 a ^ 4) + O(t^6)
                return 1 - ((1 + a2) * pt * pt) / (6 * a2);
            }
            return a * sin(pt) * sin(pt / a) / (pt * pt);
        }
    };





    // store precomputed filter weights for resampling
    struct ResamplingWeights
    {
        // max possible coeffs per output sample
        int maxCoeffs_;
        // precomputed coeffs, maxCoeffs_ of the per output sample
        vector<double> coeffs_;
        // index and delta bounds per output sample
        // todo - handle boundary conditions 
        vector<std::tuple<int, int>> bounds_;

        // enforceBoundaries allows samples to go off end, and the caller should handle them on read/write
        ResamplingWeights(
            int inSize, int outSize,
            const FilterDef& filter,
            // clamp sample indexes to range, otherwise image get/set do boundary checking
            bool enforceBoundaries,
            // perform filter size changes on downscaling to mitigate aliasing
            bool useFilterScaling
        )
        {
            // how to scale sample positions
            const double sampleScale = (double)(inSize) / outSize;
            // how to scale filter, clamp helps downsampling antialiasing
            const double filterScale =
                useFilterScaling
                ? std::max(1.0, sampleScale)
                : 1.0;
            const double invFilterScale = 1.0 / filterScale;

            // support size (length of resampling filter)
            const double filterSupport = filter.support() * filterScale;
            maxCoeffs_ = static_cast<int>(ceil(filterSupport)) * 2 + 1; // for odd filters

            coeffs_.resize(outSize * maxCoeffs_); // zeroes out entries

            // space for boundaries, each entry is a (index, delta)
            bounds_.resize(outSize);

            for (int sampleIndex = 0; sampleIndex < outSize; sampleIndex++)
            {
                // +0.5 for pixel centers
                const double center = sampleScale * (sampleIndex + 0.5);
                // Round the values
                // todo - handle boundaries here - we should let these go off the ends

                // we let these go off image ends, outer layers need to handle boundaries
                // + 0.5 for rounding (not quite correct, due to floating point issues)
                int tMin = static_cast<int>(center - filterSupport + 0.5);
                int tMax = static_cast<int>(center + filterSupport + 0.5);
                if (enforceBoundaries)
                {
                    tMin = std::max(tMin, 0);
                    tMax = std::min(tMax, inSize);
                }
                const int tDelta = tMax - tMin;
                bounds_[sampleIndex] = { tMin, tDelta };

                double* sampleKernel = &coeffs_[sampleIndex * maxCoeffs_];

                double totalWeight = 0.0; // for normalizing the kernel
                for (int x = 0; x < tDelta; x++)
                {
                    const double weight = filter.filter((x + tMin - center + 0.5) * invFilterScale);
                    sampleKernel[x] = weight;
                    totalWeight += weight;
                }
                if (totalWeight != 0.0)
                {
                    for (int x = 0; x < tDelta; x++)
                        sampleKernel[x] /= totalWeight;
                }
            }
        }
    };

    template<typename TColor, typename TImage>
    void Resample(TImage& imIn, TImage& imOut, int outWidth, int outHeight, int offset, const ResamplingWeights& coeffs, bool doVertical)
    {
        for (int j = 0; j < outHeight; j++) {
            for (int i = 0; i < outWidth; i++) {
                int x0, y0, dx, dy, tMin, tDelta;
                const double* kernel;
                if (doVertical)
                {
                    auto [s, d] = !doVertical ? coeffs.bounds_[i] : coeffs.bounds_[j];
                    tMin = s; tDelta = d;
                    x0 = i; y0 = tMin; dx = 0; dy = 1;
                    kernel = &coeffs.coeffs_[j * coeffs.maxCoeffs_];
                }
                else
                {
                    auto [s, d] = !doVertical ? coeffs.bounds_[i] : coeffs.bounds_[j];
                    tMin = s; tDelta = d;
                    x0 = tMin; y0 = j + offset; dx = 1; dy = 0;
                    kernel = &coeffs.coeffs_[i * coeffs.maxCoeffs_];
                }

                TColor ss{}; // should be all zeros
                ss = ss - ss; // ensure all zeros
                for (int t = 0; t < tDelta; t++) {
                    ss += kernel[t] * imIn.Get(x0, y0);
                    x0 += dx;
                    y0 += dy;
                }
                imOut.Set(i, j, ss);
            }
        }
    }


    ImagePtr ResampleImage(
        ImagePtr imIn,
        int w, int h,
        const FilterDef& filter,
        bool enforceBoundaries,
        bool useFilterScaling
    ) {
        int w2 = w, h2 = h;
        auto imOut = Image::Make(w2, h2);

        auto [w1, h1] = imIn->Size();

        // precompute coefficients for horizontal and vertical passes once
        ResamplingWeights horizCoeffs(
            w1, w2,
            filter,
            enforceBoundaries,
            useFilterScaling
        );
        ResamplingWeights vertCoeffs(
            h1, h2,
            filter,
            enforceBoundaries,
            useFilterScaling
        );

        // First used row in the source image
        auto [ybox_first, _] = vertCoeffs.bounds_[0];
        // Last used row in the source image
        auto [lastA, lastB] = vertCoeffs.bounds_[h2 - 1];
        int ybox_last = lastA + lastB;
        if (!enforceBoundaries)
        {
            // letting stuff out of bounds requires clamping here
            ybox_first = std::max(0, ybox_first);
            ybox_last = std::min(h1, ybox_last);
        }

        // only do subset of rows
        for (int i = 0; i < h2; i++) {
            auto [a, b] = vertCoeffs.bounds_[i];
            vertCoeffs.bounds_[i] = { a, b - ybox_first };
        }

        auto imHeight = ybox_last - ybox_first;
        auto imTemp = Image::Make(w2, imHeight);

        // horizontal pass
        Resample<Color>(*imIn, *imTemp, w2, imHeight, ybox_first, horizCoeffs, false);

        // vertical pass
        Resample<Color>(*imTemp, *imOut, w2, h2, 0, vertCoeffs, true);

        return imOut;
    }



} // namespace 