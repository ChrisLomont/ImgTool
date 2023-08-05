#pragma once
#include <vector>
#include <iostream>
#include <algorithm>
#include <numbers>


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#ifdef _MSC_VER  // needed for visual c++
#define __STDC_LIB_EXT1__ // sprintf -> sprintf_s
#endif
#include "stb_image_write.h"

#include "Timer.h"
#include "filters.h"
#include "ImageMetrics.h"
#include "Image.h"

using namespace std; // todo - remove

/*------------------- filter helpers -----------------------*/

// todo - move these, and filters, to doubles
using Run1D = vector<float>;

Run1D GetRun1D(ImagePtr img, int i0, int j0, int di, int dj, int channel)
{
	auto [w, h] = img->Size();
	Run1D run;
	int i = i0, j = j0;
	while (0 <= i && i < w && 0 <= j && j < h)
	{
		auto c = img->Get(i, j);
		run.push_back(c[channel]);
		i += di;
		j += dj;
	}
	return run;
}
void SetRun1D(ImagePtr img, int i0, int j0, int di, int dj, int channel, const Run1D& run)
{
	auto [w, h] = img->Size();
	int i = i0, j = j0, k = 0;
	while (0 <= i && i < w && 0 <= j && j < h)
	{
		auto c = img->Get(i, j);
		c[channel] = run[k++];
		img->Set(i, j, c);
		i += di;
		j += dj;
	}
}

template<typename Kernel>
vector<float> Filter1D(const vector<float>& src, int w)
{
	auto a = (int)src.size();
	Kernel k;
	if (a < w)
		return upsample(src, w, k);
	else if (a > w)
		return downsample(src, w, k);
	return src;
}

template <typename Kernel>
ImagePtr ApplyFilter(ImagePtr src, int w, int h)
{
	const auto [w1, h1] = src->Size();
	const int w2 = w, h2 = h;
	// stretch width:
	ImagePtr temp = make_shared<Image>(w2, h1);
	for (int j = 0; j < h1; ++j) // each original row
	{
		for (int channel = 0; channel < 3; ++channel) // rgb
		{
			auto run1D = GetRun1D(src, 0, j, 1, 0, channel);
			run1D = Filter1D<Kernel>(run1D, w2); // expand row
			SetRun1D(temp, 0, j, 1, 0, channel, run1D);
		}
	}
	//return temp;
	// stretch height
	ImagePtr dst = make_shared<Image>(w2, h2);
	for (int i = 0; i < w2; ++i) // each stretched column
	{
		for (int channel = 0; channel < 3; ++channel) // rgb
		{
			auto run1D = GetRun1D(temp, i, 0, 0, 1, channel);
			run1D = Filter1D<Kernel>(run1D, h2);
			SetRun1D(dst, i, 0, 0, 1, channel, run1D);
		}
	}
	return dst;
}


/* -------------------- Resizers ------------------------------*/
ImagePtr ResizeNN(ImagePtr src, int w, int h)
{
	auto [w1, h1] = src->Size();
	ImagePtr dst = make_shared<Image>(w, h, 4);
	// iterate over dst
	for (int j = 0; j < h; ++j)
	{
		double jp = j + 0.5; // pixel centers
		jp *= h1;
		jp /= h; // now in src coords
		// src coord
		int sj = (int)clamp(round(jp - 0.5), 0.0, h1 - 1.0);

		for (int i = 0; i < w; ++i)
		{
			double ip = i + 0.5; // pixel centers
			ip *= w1;
			ip /= w; // now in src coords
			// src coord
			int si = (int)clamp(round(ip - 0.5), 0.0, w1 - 1.0);
			dst->Set(i, j, src->Get(si, sj));
		}
	}
	return dst;
}
ImagePtr ResizeBilinear(ImagePtr src, int w, int h)
{
	// return ApplyFilter<Box>(src, w, h);
	return ApplyFilter<Hat>(src, w, h);
}
ImagePtr ResizeBicubic(ImagePtr src, int w, int h)
{
	// also known as Keys
	return ApplyFilter<CatmullRom>(src, w, h);
}


template<size_t N1>
struct Lanczos : KernelBase<2 * N1> {
	float operator()(float x) const override {
		// Lanczos size -a <= x <= a
		//constexpr float pi = 3.14159265358979323846f;
		constexpr float a = N1;// 0.5f * support();
		if (x <= -a || a <= x) return 0;
		if (fabs(x) < 1e-5) // todo - compute cutoff correctly
		{
			// taylor series from mathematica Series[Sinc[\[Pi] x] Sinc[\[Pi] x / a], {x, 0, 5}]
			//1 - (((1 + a^2) \[Pi]^2) x^2)/(6 a ^ 2) + ((3 + 10 a ^ 2 + 3 a ^ 4) \[Pi] ^ 4 x ^ 4) / (360 a ^ 4)
				// todo - precompute constants, then simply mult here
			constexpr float a2 = a * a;
			constexpr float pi2 = (float)(numbers::pi * numbers::pi);
			constexpr float a4 = a2 * a2;
			constexpr float pi4 = pi2 * pi2;
			constexpr float c1 = (1 + a2) * pi2 / (6 * a2);
			constexpr float c2 = (3 + 10 * a2 + 3 * a4) * pi4 / (360 * a4);
			float x2 = x * x;
			float x4 = x2 * x2;
			return 1 - c1 * x2 + c2 * x4;
		}
		auto px = numbers::pi * x;
		return a * sin(px) * sin(px / a) / (px * px);
	}
	void accumulate_buffer(float fu, float u) override {
		throw exception();
	}
	float sample_buffer(float u) const override {
		throw exception();
	}
	void digital_filter(vector<float>&) const override { /* todo - anything here? */ }
	string name() const override { return string("Lanczos") + to_string(N1); }

	// How much does the kernel integrate to?
	// obtained from Mathematica
	float integral() const override {
		static float intg[]{
			1.009789840618732f, // size 2
			0.9970553459543397f, // size 3
			1.001256108101727f, // size 4
			0.9993526254518625f,
			1.000376158906594f
		};
		if (N1 < 2 || 6 < N1) throw exception();
		return intg[N1 - 2];

	}
};
struct  Lanczos2 : Lanczos<2> {};
struct  Lanczos3 : Lanczos<3> {};
struct  Lanczos4 : Lanczos<4> {};
struct  Lanczos5 : Lanczos<5> {};
struct  Lanczos6 : Lanczos<6> {};

double LanczosFunc(double v, int a)
{
	v = abs(v);
	if (v >= a) return 0;
	if (v == 0.0) return 1;
	double px = numbers::pi * v;
	// todo - for near 0, use taylor series, compute a few and put here
	return a * sin(px) * sin(px / a) / (px * px);
}

ImagePtr ResizeLanczosRadial(ImagePtr src, int w, int h, double a)
{ // radial lanczos
	auto [w1, h1] = src->Size();
	int ws = w1, hs = h1; // clang error!
	auto dst = make_shared<Image>(w,h);
	int min = (int)floor(a), max = (int)ceil(a);
	auto r2 = a * a;
	dst->Apply(
		[&](int i, int j)
		{
			// center pixel into src location
			double cx = ((i + 0.5) * ws) / w;
			double cy = ((j + 0.5) * hs) / h;
			double sum = 0;
			Color c(0,0,0,1);
			for (int sj = (int)floor(cy-a); sj <= (int)ceil(cy+a); ++sj)
				for (int si = floor(cx - a); si <= ceil(cx + a); ++si)
				{
					//	if (!src->Legal(si, sj)) continue; 
					auto dx = cx - (si + 0.5);
					auto dy = cy - (sj + 0.5);
					auto d2 = dx * dx + dy * dy;
					if (d2 <= r2)
					{
						auto wt = LanczosFunc(sqrt(d2),(int)a);
						sum += wt;
						c += wt * src->Get(si, sj);						
					}					
				}
			
			c /= (sum>0)?sum:1;
			c.a = 1.0; // todo?
			return c;
		}
	);

	return dst;	
}



ImagePtr ResizeLanczos(ImagePtr src, int w, int h, double a)
{
	// todo; - make using kernel?
	const auto [w1, h1] = src->Size();
	const int w2 = w, h2 = h;
	// stretch width:
	ImagePtr temp = make_shared<Image>(w2, h1);
	for (int j = 0; j < h1; ++j) // each temp row
		for (int i = 0; i < w2; ++i) // each temp column
		{
			// compute source index of pixel center:
			double px = (i + 0.5) * w1 / w2;
			int i0 = (int)(floor(px - a - 1));
			int i1 = i0 + 2 * a;
			Color c(0, 0, 0, 0);
			double total = 0.0;
			for (int ii = i0; ii <= i1; ++ii)
			{
				double dx = (ii + 0.5) - px; // dist in source space
				double weight = LanczosFunc(dx, a);
				total += weight;
				c += weight * src->Get(ii, j);
			}
			c = 1.0/total * c; // todo - div by 0 ?
			c.a = 1; // todo - blend?
			temp->Set(i, j, c);
		}

	// stretch height
	ImagePtr dst = make_shared<Image>(w2, h2);

	for (int j = 0; j < h2; ++j) // each dest row
		for (int i = 0; i < w2; ++i) // each dest column
		{
			// compute source index of pixel center:
			double py = (j + 0.5) * h1 / h2;
			int j0 = (int)(floor(py - a - 1));
			int j1 = j0 + 2 * a;
			Color c(0, 0, 0, 0);
			double total = 0.0;
			for (int jj = j0; jj <= j1; ++jj)
			{
				double dy = (jj + 0.5) - py; // dist in source space
				double weight = LanczosFunc(dy, a);
				total += weight;
				c += weight * temp->Get(i, jj);
			}
			c = 1.0 / total * c; // todo - div by 0?
			c.a = 1; // todo - blend?
			dst->Set(i, j, c);
		}
	return dst;
}

// alpha blend
Color Blend(const Color & over, const Color & under)
{
	auto alphaOver = over.a;
	auto alphaUnder = under.a;
	auto alpha = alphaOver + alphaUnder * (1 - alphaOver);
	//c = (cOver*aOver + cUnder*bUnder*(1-aover))/alpha;
	auto c = (alphaOver * over + alphaUnder * (1 - alphaOver) * under);
	if (alpha != 0)
		c = (1.0 / alpha) * c;
	else
		c.r = c.g = c.b = c.a = 0;
	c.a = alpha;
	return c;
}

void Blit(ImagePtr dst, int dx, int dy, ImagePtr src, int x1, int y1, int w, int h)
{
	for (int j = 0; j < h; ++j)
		for (int i = 0; i < w; ++i)
		{
			auto srcColor = src->Get(i + x1, j + y1);
			auto dstColor = dst->Get(i + dx, j + dy);
			auto color = Blend(srcColor, dstColor);
			dst->Set(i + dx, j + dy, color);
		}
}

/* -------------------- Metrics ------------------------------*/
double MetricMSE(ImagePtr src, ImagePtr dst)
{
	auto [w, h] = dst->Size();
	return Lomont::Graphics::ImageMetrics::MSE(
		w, h,
		[&](int i, int j) {return src->Get(i, j).r; },
		[&](int i, int j) {return dst->Get(i, j).r; }
	);
}
double MetricPSNR(ImagePtr src, ImagePtr dst)
{
	auto [w, h] = dst->Size();
	return Lomont::Graphics::ImageMetrics::PSNR(
		w, h,
		[&](int i, int j) {return src->Get(i, j).r; },
		[&](int i, int j) {return dst->Get(i, j).r; }
	);
}
double MetricSSIM(ImagePtr src, ImagePtr dst)
{
	auto [w, h] = dst->Size();
	return Lomont::Graphics::ImageMetrics::SSIM(
		w, h,
		[&](int i, int j) {return src->Get(i, j)[0]; },
		[&](int i, int j) {return dst->Get(i, j)[0]; }
	);
}



/* -------------------- Rotation -----------------------------*/

typedef std::function<Color(const ImagePtr&, double, double)> InterpFunc;

// loop over dest pixels, backproject, interp
ImagePtr RotateDest(const ImagePtr& src, double angleInRadians, bool expand, InterpFunc Interp)
{
	auto [ws, hs] = src->Size();

	// treat pixel centers as 0.5, 1.5, 2.5, ... w-0.5

	// get center of src image
	double csx = ws / 2.0;
	double csy = hs / 2.0;

	int radius, sideLen;
	double cdx, cdy;
	if (expand)
	{

		// compute center of dest image
		// want integer or half integer the same as src image to align pixels better
		// embed in circle inside square image
		radius = ceil(sqrtf(csx * csx + csy * csy)); // at least this
		cdx = radius + (ws & 1) / 2.0f;
		cdy = radius + (hs & 1) / 2.0f;
		sideLen = ceil(2 * max(cdx, cdy)); // big enough
	}
	else
	{
		sideLen = max(ws, hs);
		cdx = cdy = sideLen / 2.0f;
	}

	auto dst = Image::Make(sideLen, sideLen);

	// angle backwards to invert projection
	double cc = cos(angleInRadians);
	double ss = sin(angleInRadians);

	// iterate over dest pixels
	for (int i = 0; i < sideLen; ++i)
		for (int j = 0; j < sideLen; ++j)
		{
			// distance from dest pixel center to center of dest image
			double dx = (i + 0.5) - cdx;
			double dy = (j + 0.5) - cdy;

			//if (i2*i2 + j2*j2 >= r * r)
			//	continue; // not in frame


			// position in source image (rotate dx,dy)
			double sx = cc * dx - ss * dy + csx;
			double sy = ss * dx + cc * dy + csy;

			// now do interpolation:
			Color c = Interp(src, sx, sy);
			c.Clamp();
			dst->Set(i, j, c);
		}

	return dst;
}
Color Interp(const Color& c1, const Color& c2, double interp)
{
	return (1 - interp) * c1 + interp * c2;
}

// interpolate the pixel location
Color InterpNN(const ImagePtr& src, double x, double y)
{
	int si = round(x);
	int sj = round(y);
	return src->Get( si, sj);
}
// interpolate the pixel location
Color InterpBilinear(const ImagePtr& src, double x, double y)
{
	// todo - this not quite bilin? check carefully
	int si = floor(x);
	int sj = floor(y);
	double fi = x - si;
	double fj = y - sj;

	auto c00 = src->Get(si, sj);
	auto c10 = src->Get( si + 1, sj);

	auto c01 = src->Get(si, sj + 1);
	auto c11 = src->Get(si + 1, sj + 1);

	auto color = Interp(
		Interp(c00, c10, fi),
		Interp(c01, c11, fi),
		fj);
	color.a = 1;
	return color;
}


// todo - allow other types, use mitchell B,C parameterization?
// for t in 0-1
// interpolate the pixel location, default bad color outside
// uses Keys bicubic, a = -0.5, most common (?)
Color InterpBicubic(const ImagePtr& src, double x, double y)
{
	// https://en.wikipedia.org/wiki/Bicubic_interpolation

	// todo - this not quite bicub? check carefully
	int si = floor(x);
	int sj = floor(y);
	double tx = x - si; // in 0-1
	double ty = y - sj;

	// interp
	auto p = [](double t, Color c0, Color c1, Color c2, Color c3) {
		/* // todo - make matrix forms?
		1/2 * [1 t t^2 t^3]*M*[c0 c1 c2 c3]^T

			|  0  2  0  0 |
		M = | -1  0  1  0 |
			|  2 -5  4 -1 |
			| -1  3 -3  1 |
		*/

		// mathematica, without the /2 part
		// c2(t + 4 t ^ 2 - 3 t ^ 3) + c0(-t + 2 t ^ 2 - t ^ 3) + c3(-t ^ 2 + t ^ 3) + c1(2 - 5 t ^ 2 + 3 t ^ 3)
		assert(0 <= t && t < 1);
		auto t2 = t * t;
		auto t3 = t * t2;
		auto c =
			c0 * (-t + 2 * t2 - t3) +
			c1 * (2 - 5 * t2 + 3 * t3) +
			c2 * (t + 4 * t2 - 3 * t3) +
			c3 * (-t2 + t3);
		c = c * 0.5f;
		c.a = 1.0; // todo - right?
		c.Clamp();
		return c;
	};
	// sample
	auto f = [&](int di, int dj) {
		return src->Get( si + di, sj + dj);
	};

	// bn = b_{-1}
	auto bn = p(tx, f(-1, -1), f(+0, -1), f(+1, -1), f(2, -1));
	auto b0 = p(tx, f(-1, +0), f(+0, +0), f(+1, +0), f(2, +0));
	auto b1 = p(tx, f(-1, +1), f(+0, +1), f(+1, +1), f(2, +1));
	auto b2 = p(tx, f(-1, +2), f(+0, +2), f(+1, +2), f(2, +2));
	return p(ty, bn, b0, b1, b2);
}



/* -------------------- Shift --------------------------------*/

ImagePtr ShiftImage2(ImagePtr img, double dx, double dy, InterpFunc Interp)
{
	auto [w, h] = img->Size();
	auto dst = make_shared<Image>(w,h);
	int w1 = w, h1 = h; // avoid clang error

	dst->Apply([&](int i, int j) {
		return img->Get(i+dx,j+dy);
		});
	return dst;
}

