#pragma once
#include <vector>
#include <iostream>
#include <algorithm>
#include <numbers>




#include "filters.h"
#include "ImageMetrics.h"
#include "Image.h"

using namespace std; // todo - remove

/*------------------- filter helpers -----------------------*/

// todo - move these, and filters, to doubles
using Run1D = vector<float>;

inline Run1D GetRun1D(ImagePtr img, int i0, int j0, int di, int dj, int channel)
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
	const auto a = static_cast<int>(src.size());
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
	const ImagePtr temp = Image::Make(w2, h1);
	for (int j = 0; j < h1; ++j) // each original row
	{
		for (int channel = 0; channel < 4; ++channel) // rgba
		{
			auto run1D = GetRun1D(src, 0, j, 1, 0, channel);
			run1D = Filter1D<Kernel>(run1D, w2); // expand row
			SetRun1D(temp, 0, j, 1, 0, channel, run1D);
		}
	}
	//return temp;
	// stretch height
	ImagePtr dst = Image::Make(w2, h2);
	for (int i = 0; i < w2; ++i) // each stretched column
	{
		for (int channel = 0; channel < 4; ++channel) // rgba
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
	ImagePtr dst = Image::Make(w, h);
	// iterate over dst
	for (int j = 0; j < h; ++j)
	{
		double jp = j + 0.5; // pixel centers
		jp *= h1;
		jp /= h; // now in src coords
		// src coord
		const int sj = static_cast<int>(clamp(round(jp - 0.5), 0.0, h1 - 1.0));

		for (int i = 0; i < w; ++i)
		{
			double ip = i + 0.5; // pixel centers
			ip *= w1;
			ip /= w; // now in src coords
			// src coord
			const int si = static_cast<int>(clamp(round(ip - 0.5), 0.0, w1 - 1.0));
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
			constexpr float pi2 = static_cast<float>(numbers::pi * numbers::pi);
			constexpr float a4 = a2 * a2;
			constexpr float pi4 = pi2 * pi2;
			constexpr float c1 = (1 + a2) * pi2 / (6 * a2);
			constexpr float c2 = (3 + 10 * a2 + 3 * a4) * pi4 / (360 * a4);
			const float x2 = x * x;
			const float x4 = x2 * x2;
			return 1 - c1 * x2 + c2 * x4;
		}
		const auto px = numbers::pi * x;
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
	if (v < 1e-20) return 1; // todo - analyze cutoff
	const double px = numbers::pi * v;
	// todo - for near 0, use taylor series, compute a few and put here
	return a * sin(px) * sin(px / a) / (px * px);
}

ImagePtr ResizeLanczosRadial(ImagePtr src, int w, int h, double a)
{ // radial lanczos
	auto [w1, h1] = src->Size();
	const int ws = w1, hs = h1; // clang error!
	auto dst = Image::Make(w,h);
	int min = static_cast<int>(floor(a)), max = static_cast<int>(ceil(a));
	const auto r2 = a * a;
	dst->Apply(
		[&](int i, int j)
		{
			// center pixel into src location
			const double cx = ((i + 0.5) * ws) / w;
			const double cy = ((j + 0.5) * hs) / h;
			double sum = 0;
			Color c(0,0,0,1);
			for (int sj = static_cast<int>(floor(cy - a)); sj <= static_cast<int>(ceil(cy + a)); ++sj)
				for (int si = floor(cx - a); si <= ceil(cx + a); ++si)
				{
					//	if (!src->Legal(si, sj)) continue; 
					const auto dx = cx - (si + 0.5);
					const auto dy = cy - (sj + 0.5);
					const auto d2 = dx * dx + dy * dy;
					if (d2 <= r2)
					{
						const auto wt = LanczosFunc(sqrt(d2),static_cast<int>(a));
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



ImagePtr ResizeLanczos(ImagePtr src, int w, int h, int a)
{
	// todo; - make using kernel?
	const auto [w1, h1] = src->Size();
	const int w2 = w, h2 = h;
	// stretch width:
	const ImagePtr temp = Image::Make(w2, h1);
	for (int j = 0; j < h1; ++j) // each temp row
		for (int i = 0; i < w2; ++i) // each temp column
		{
			// compute source index of pixel center:
			const double px = (i + 0.5) * w1 / w2;
			const int i0 = static_cast<int>(floor(px - a - 2));
			const int i1 = i0 + 2 * a+1;
			Color c(0, 0, 0, 0);
			double total = 0.0;
			for (int ii = i0; ii <= i1; ++ii)
			{
				const double dx = (ii + 0.5) - px; // dist in source space
				const double weight = LanczosFunc(dx, a);
				total += weight;
				c += weight * src->Get(ii, j);
			}
			c /= total; // todo - div by 0 ?			
			c.Clamp();
			temp->Set(i, j, c);
		}

	// stretch height
	ImagePtr dst = Image::Make(w2, h2);

	for (int j = 0; j < h2; ++j) // each dest row
		for (int i = 0; i < w2; ++i) // each dest column
		{
			// compute source index of pixel center:
			const double py = (j + 0.5) * h1 / h2;
			const int j0 = static_cast<int>(floor(py - a - 2));
			const int j1 = j0 + 2 * a+1;
			Color c(0, 0, 0, 0);
			double total = 0.0;
			for (int jj = j0; jj <= j1; ++jj)
			{
				const double dy = (jj + 0.5) - py; // dist in source space
				const double weight = LanczosFunc(dy, a);
				total += weight;
				c += weight * temp->Get(i, jj);
			}
			c /= total; // todo - div by 0?			
			c.Clamp();
			dst->Set(i, j, c);
		}
	return dst;
}

// alpha blend
// assumes pre-multiplied alpha
// does Porter-Duff 1984 OVER operation
Color Blend(const Color& under, const Color & over)
{
	// in pre-multiplied alpha, OVER is a 4-vector operation:
	// https://en.wikipedia.org/wiki/Alpha_compositing
	return over + under * (1 - over.a);
}

// pre-multiply alpha or reverse it
inline void AlphaCorrect(ImagePtr img, bool premultiplyAlpha)
{
	if (premultiplyAlpha)
		img->Apply([](const Color& color) { return color.ToPremultipliedAlpha(); });
	else
		img->Apply([](const Color& color) { return color.FromPremultipliedAlpha(); });
}


inline ImagePtr Blit(
	const ImagePtr& underImage, int dx, int dy, ImagePtr overImage, 
	const int x1, int y1, int w, int h, 
	bool alphaBlend)
{
	auto result = Image::Make(underImage);
	Color color;
	for (int j = 0; j < h; ++j)
		for (int i = 0; i < w; ++i)
		{
			auto overColor = overImage->Get(i + x1, j + y1);
			if (!alphaBlend)
			{
				color = overColor;
			}
			else 
			{
				auto underColor = underImage->Get(i + dx, j + dy);
				color = Blend(underColor, overColor);
			}
			result->Set(i + dx, j + dy, color);
		}
	return result;
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

/* crop */
ImagePtr CropImage(const ImagePtr& src, int x1, int y1, int x2, int y2)
{
	if (x2 < x1) swap(x1, x2);
	if (y2 < y1) swap(y1, y2);
	auto [w, h] = src->Size();
	int w2 = x2 - x1 + 1, h2 = y2 - y1 + 1;
	const auto dst = Image::Make(w2, h2);
	for (int j = 0; j < h2; ++j)
		for (int i = 0; i < w2; ++i)
			dst->Set(i, j, src->Get(i + x1, j + y1));
	return dst;
}


/* -------------------- Rotation -----------------------------*/

typedef std::function<Color(const ImagePtr&, double, double)> InterpFunc;

// loop over dest pixels, backproject, interp
ImagePtr RotateDest(const ImagePtr& src, double angleInRadians, bool expand, InterpFunc Interp)
{
	auto [ws, hs] = src->Size();

	// treat pixel centers as 0.5, 1.5, 2.5, ... w-0.5

	// get center of src image
	const double csx = ws / 2.0;
	const double csy = hs / 2.0;

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
	const double cc = cos(angleInRadians);
	const double ss = sin(angleInRadians);

	// iterate over dest pixels
	for (int i = 0; i < sideLen; ++i)
		for (int j = 0; j < sideLen; ++j)
		{
			// distance from dest pixel center to center of dest image
			const double dx = (i + 0.5) - cdx;
			const double dy = (j + 0.5) - cdy;

			//if (i2*i2 + j2*j2 >= r * r)
			//	continue; // not in frame


			// position in source image (rotate dx,dy)
			const double sx = cc * dx - ss * dy + csx;
			const double sy = ss * dx + cc * dy + csy;

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
	const int si = round(x);
	const int sj = round(y);
	return src->Get( si, sj);
}
// interpolate the pixel location
Color InterpBilinear(const ImagePtr& src, double x, double y)
{
	// todo - this not quite bilin? check carefully
	const int si = floor(x);
	const int sj = floor(y);
	const double fi = x - si;
	const double fj = y - sj;

	const auto c00 = src->Get(si, sj);
	const auto c10 = src->Get( si + 1, sj);

	const auto c01 = src->Get(si, sj + 1);
	const auto c11 = src->Get(si + 1, sj + 1);

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
	const int si = floor(x);
	const int sj = floor(y);
	const double tx = x - si; // in 0-1
	const double ty = y - sj;

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
		const auto t2 = t * t;
		const auto t3 = t * t2;
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
	const auto bn = p(tx, f(-1, -1), f(+0, -1), f(+1, -1), f(2, -1));
	const auto b0 = p(tx, f(-1, +0), f(+0, +0), f(+1, +0), f(2, +0));
	const auto b1 = p(tx, f(-1, +1), f(+0, +1), f(+1, +1), f(2, +1));
	const auto b2 = p(tx, f(-1, +2), f(+0, +2), f(+1, +2), f(2, +2));
	return p(ty, bn, b0, b1, b2);
}


ImagePtr ShearX(const ImagePtr& img, float shear)
{
	auto [w1, h1] = img->Size();
	int shift = round(h1 * shear);
	int w2 = w1 + abs(shift) + 1;
	if (shift < 0) shift = 0;
	auto dst = Image::Make(w2, h1);

	for (int row = 0; row < h1; ++row)
	{
		auto r0 = img->Row(row);
		auto r1 = dst->Row(row);
		int offset = shift - shear * row; // todo - this important :)
		for (int col = 0; col < w1; ++col)
		{
			assert(0 <= col && col < r0.num);
			r1[col + offset] = r0[col];
		}
	}
	return dst;
}
// rotate by given number of 90 degree rotations
ImagePtr RotateOrtho(const ImagePtr& img, int rots90, bool expand = false)
{
	auto [w, h] = img->Size();

	rots90 %= 4; // -3 to 3
	rots90 = (rots90 + 4) % 4; // 0-3

	bool odd = (rots90 & 1) != 0;

	auto dst = Image::Make(odd ? h : w, odd ? w : h);

	dst->Apply([&](int i, int j)
		{
			Color c(0, 0, 0, 0);
			if (rots90 == 0) c = img->Get(i, j);
			if (rots90 == 1) c = img->Get(w - 1 - j, i);
			if (rots90 == 2) c = img->Get(w - 1 - i, h - 1 - j);
			if (rots90 == 3) c = img->Get(j, h - 1 - i);
			return c;
		}
	);
	return dst;
}
ImagePtr ShearY(const ImagePtr& img, float shear)
{
	// todo- make faster
	auto rot1 = RotateOrtho(img, -1);
	auto sh = ShearX(rot1, shear);
	return RotateOrtho(sh, +1);
}

extern ImagePtr CropToBounds(const ImagePtr& img, const Color& c );

ImagePtr Rotate3Shear(const ImagePtr& img, float angleInRadians, bool expand)
{
	// jupyter nb at https://github.com/str4w/ImageRotation_PWL/blob/master/UnserRotation.ipynb
	// shear in x by -tan(angle/2) (explodes!)
	// shear in y by sin(angle)
	// shear in x by -tan(angle/2) (again! explodes!)
	// equivalent angle in [0,2pi)
	// also implement FULLY REVERSIBLE IMAGE ROTATION BY 1-D FILTERING - 2008 - better filter, works nicely

	auto a1 = fmodf(angleInRadians, numbers::pi * 2); // -2pi to +2pi
	while (a1 < 0) a1 += 2 * numbers::pi; // note 'if' may not be enough? (todo - check carefully)
	assert(a1 >= 0);
	auto q1 = fmodf(a1, numbers::pi / 2); // excess into quadrant
	assert(q1 >= 0);
	int rot90 = roundf((a1 - q1) / (numbers::pi / 2)); // num of rots 90

	// final fixup: instead of 0-pi/2 rotation, make it [-pi/4,pi/4]
	if (q1 > numbers::pi / 4)
	{
		rot90++;
		q1 -= numbers::pi / 2;
	}
	const ImagePtr src = (rot90 > 0) ? RotateOrtho(img, rot90) : img;

	auto shearX = -tanf(q1 / 2);
	auto shearY = sinf(q1);
	auto s1 = ShearX(src, shearX);
	auto s2 = ShearY(s1, -shearY);
	auto s3 = ShearX(s2, shearX);
	if (expand)
		return s3;	
	return CropToBounds(s3, Color(0,0,0)); // todo - image was too large, clip more intelligently
}


/*** DCT Rotation ***/

void DCT(const vslice<Color>& src, vslice<Color>& dst, size_t n)
{ // DCT-II on wikipedia
	float sc = numbers::pi / n;
	for (int k = 0; k < n; ++k)
	{
		dst[k] = Color(0, 0, 0);
		for (int j = 0; j < n; ++j)
			dst[k] += src[j] * cosf(sc * (j + 0.5f) * k);
	}
}
// eval IDCT at given value t
Color IDCT(const vslice<Color>& src, float t, size_t n)
{ // DCT-II on wikipedia
	float invn = 1.0f / n;
	float sc = numbers::pi * invn;
	auto c = src[0];
	for (int j = 1; j < n; ++j)
		c += src[j] * 2 * cosf(sc * j * (t + 0.5));
	return c * invn;
}


// eval at s,t
Color EvalIDCT(const ImagePtr& src, float s, float t, int n)
{ // DCT-III
	vector<Color> tmp(n);
	for (int row = 0; row < n; ++row)
	{
		tmp[row] = IDCT(src->Row(row), s, n);
	}
	auto col = vslice(tmp);
	auto color = IDCT(col, t, n);
	// clean color
	color.a = 1.0;
	color.Clamp();
	return color;
}


Color InterpDCT(const ImagePtr& src, double x, double y)
{   
	// todo - check carefully
	int si = (int)floor(x);
	int sj = (int)floor(y);
	double tx = x - si; // in 0-1
	double ty = y - sj;

	// 9x9 size
	const int n = 9;
	const ImagePtr f = Image::Make(n, n);
	ImagePtr g = Image::Make(n, n);
	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)
			f->Set(i, j, src->Get(si + i - n / 2, sj + j - n / 2));
	for (int row = 0; row < n; ++row)
	{
		auto gr = g->Row(row); // cannot ref temp var
		DCT(f->Row(row), gr, n);
	}
	for (int col = 0; col < n; ++col)
	{
		auto fc = f->Col(col); // cannot ref temp var
		DCT(g->Col(col), fc, n);
	}

	// eval the inverse at tx,ty
	return EvalIDCT(f, tx + n / 2, ty + n / 2, n);
}



/* -------------------- Shift --------------------------------*/

ImagePtr ShiftImage2(ImagePtr img, double dx, double dy, InterpFunc Interp)
{
	auto [w, h] = img->Size();
	auto dst = Image::Make(w,h);
	int w1 = w, h1 = h; // avoid clang error

	dst->Apply([&](int i, int j) {
		return img->Get(i+dx,j+dy);
		});
	return dst;
}

