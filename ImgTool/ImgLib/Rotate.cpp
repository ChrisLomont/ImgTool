// lomont old rotate code
#if 0
// Paeth style rotation: 3 shears
// see https://gautamnagrawal.medium.com/rotating-image-by-any-angle-shear-transformation-using-only-numpy-d28d16eb5076
// https://datagenetics.com/blog/august32013/index.html
// http://archive.gamedev.net/archive/reference/articles/article811.html
#include <cmath>
#include <cassert>
#include "../Image/Image.h"
#include "Rotate.h"
#include "../Image/Draw.h"
#include "../Base/Utils.h"
using namespace std;

// todo - make samplers that wrap, clamp, reflect, etc.
Color Sample(const ImagePtr& src, int i, int j)
{
	auto [w, h] = src->Size();
	i = clamp(i, 0, w - 1);
	j = clamp(j, 0, h - 1);
	//if (0 <= i && 0 <= j && i < w && j < h)
	return src->Get(i, j);
	return Color(1, 0, 1); // todo - reflect?
}

// todo - template
Color Interp(const Color& c1, const Color& c2, float interp)
{
	return c1 * (1 - interp) + c2 * interp;
}

// interpolate the pixel location
Color InterpNN(const ImagePtr& src, int w, int h, float x, float y)
{
	int si = round(x);
	int sj = round(y);
	return Sample(src, si, sj);
}

// interpolate the pixel location
Color InterpBilin(const ImagePtr& src, int w, int h, float x, float y)
{
	// todo - this not quite bilin? check carefully
	int si = floor(x);
	int sj = floor(y);
	float fi = x - si;
	float fj = y - sj;

	auto c00 = Sample(src, si, sj);
	auto c10 = Sample(src, si + 1, sj);

	auto c01 = Sample(src, si, sj + 1);
	auto c11 = Sample(src, si + 1, sj + 1);

	auto color = Interp(
		Interp(c00, c10, fi),
		Interp(c01, c11, fi),
		fj);
	color.a = 1;
	return color;
}

// interpolate the pixel location
Color InterpShiftedBilin(const ImagePtr& src, int w, int h, float x, float y, int style)
{
	// implements "Linear Interpolation Revisited", Blu, Unser, 2004 with shift they derived of (1-sqrtf(3)/3)/2; // from paper, ~0.21
	// also can do "Shifting Interpolation Kernel Toward Orthogonal Projection" 2017, Sadeghi, Yu, Wang, with 1/8 = 0.125


	auto tau = style == 0
		? (1 - sqrtf(3) / 3) / 2 // from paper, ~0.21
		: 0.125f; // other paper

	// todo - this not quite bilin? check carefully
	int si = floor(x);
	int sj = floor(y);
	float fi = x - si + tau;
	float fj = y - sj + tau;

	auto c00 = Sample(src, si, sj);
	auto c10 = Sample(src, si + 1, sj);

	auto c01 = Sample(src, si, sj + 1);
	auto c11 = Sample(src, si + 1, sj + 1);

	auto color = Interp(
		Interp(c00, c10, fi),
		Interp(c01, c11, fi),
		fj);
	color.a = 1;
	return color;
}

// interpolate the pixel location
Color InterpShiftedBilinA(const ImagePtr& src, int w, int h, float x, float y)
{
	return InterpShiftedBilin(src, w, h, x, y, 0);
}
// interpolate the pixel location
Color InterpShiftedBilinB(const ImagePtr& src, int w, int h, float x, float y)
{
	return InterpShiftedBilin(src, w, h, x, y, 1);
}


// for t in 0-1

// interpolate the pixel location, default bad color outside
// uses Keys bicubic, a = -0.5, most common (?)
Color InterpBicub(const ImagePtr& src, int w, int h, float x, float y)
{
	// https://en.wikipedia.org/wiki/Bicubic_interpolation

	// todo - this not quite bicub? check carefully
	int si = floor(x);
	int sj = floor(y);
	float tx = x - si; // in 0-1
	float ty = y - sj;

	// interp
	auto p = [](float t, Color c0, Color c1, Color c2, Color c3) {
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
		return Sample(src, si + di, sj + dj);
	};

	// bn = b_{-1}
	auto bn = p(tx, f(-1, -1), f(+0, -1), f(+1, -1), f(2, -1));
	auto b0 = p(tx, f(-1, +0), f(+0, +0), f(+1, +0), f(2, +0));
	auto b1 = p(tx, f(-1, +1), f(+0, +1), f(+1, +1), f(2, +1));
	auto b2 = p(tx, f(-1, +2), f(+0, +2), f(+1, +2), f(2, +2));
	return p(ty, bn, b0, b1, b2);
}
#if 1
void DCT(const vslice<Color>& src, vslice<Color>& dst, size_t n)
{ // DCT-II on wikipedia
	float sc = pi / n;
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
	float sc = pi * invn;
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
#endif

// todo 
Color InterpDCT(const ImagePtr& src, int w, int h, float x, float y)
{ // fast DCT http://fourier.eng.hmc.edu/e161/lectures/dct/node2.html
	// see https://stackoverflow.com/questions/22768869/i-am-looking-for-a-simple-algorithm-for-fast-dct-and-idct-of-matrix-nxm
	// implement fast DCT, often called FDT (see wikipedia DCT)
	// https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms

	// todo - check carefully
	int si = floor(x);
	int sj = floor(y);
	float tx = x - si; // in 0-1
	float ty = y - sj;

	// 9x9 size
	const int n = 9;
#if 1
	ImagePtr f = Image::Make(n, n);
	ImagePtr g = Image::Make(n, n);
	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)
			f->Set(i, j, Sample(src, si + i - n / 2, sj + j - n / 2));
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
#else
	static vector<float> A; // TODO - this not thread safe
	auto c = [](float u) {return (u == 0) ? sqrtf(1.0f / n) : sqrtf(2.0f / n); };
	if (A.size() == 0)
	{ // compute A matrix
		A.resize(n * n);
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i)
				A[i + j * n] = c(i) * cosf((0.5f + j) * pi * i / n);
	}


	vector<Color> F(n * n);

	// DCT of nxn nbhd computed as F = A f A^T where f is nbhd
	// 1. F <- A f
	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)
		{
			Color sum(0, 0, 0, 1); // start
			for (int k = 0; k < n; ++k)
				sum += Sample(src, k - n / 2 + si, j - n / 2 + sj) * A[i + k * n];
			F[i + j * n] = sum;
		}
	// 2. F <- F * A^T
	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)
		{
			Color sum(0, 0, 0, 1); // start
			for (int k = 0; k < n; ++k)
				sum += F[i, k] * A[j + k * n];
			F[i + j * n] = sum;
		}

	// invert DCT to compute at tx,ty
	tx += n / 2; // recenter
	ty += n / 2;
	Color fij(0, 0, 0, 1);
	for (int u = 0; u < n; ++u)
		for (int v = 0; v < n; ++v)
		{
			fij += F[u + v * n] * c(u) * c(v) * cosf((0.5f + tx) * pi * u / n) * cosf((0.5f + ty) * pi * v / n);
		}
	fij.a = 1.0f;
	fij.Clamp();
	return fij;
#endif

}


typedef std::function<Color(const ImagePtr&, int, int, float, float)> InterpFunc;

// loop over dest pixels, backproject, interp
ImagePtr RotateDest(const ImagePtr& src, float angleInRadians, bool expand, InterpFunc Interp)
{
	auto [ws, hs] = src->Size();

	// treat pixel centers as 0.5, 1.5, 2.5, ... w-0.5

	// get center of src image
	float csx = ws / 2.0f;
	float csy = hs / 2.0f;

	int radius, sideLen;
	float cdx, cdy;
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
	float cc = cosf(angleInRadians);
	float ss = sinf(angleInRadians);

	// iterate over dest pixels
	for (int i = 0; i < sideLen; ++i)
		for (int j = 0; j < sideLen; ++j)
		{
			// distance from dest pixel center to center of dest image
			float dx = (i + 0.5f) - cdx;
			float dy = (j + 0.5f) - cdy;

			//if (i2*i2 + j2*j2 >= r * r)
			//	continue; // not in frame


			// position in source image (rotate dx,dy)
			float sx = cc * dx - ss * dy + csx;
			float sy = ss * dx + cc * dy + csy;

			// now do interpolation:
			Color c = Interp(src, ws, hs, sx, sy);
			c.Clamp();
			dst->Set(i, j, c);
		}

	return dst;
}

// rotate by given number of 90 degree rotations
ImagePtr RotateOrtho(const ImagePtr& img, int rots90, bool expand)
{
	auto [w, h] = img->Size();

	rots90 %= 4; // -3 to 3
	rots90 = (rots90 + 4) % 4; // 0-3

	bool odd = (rots90 & 1) != 0;

	auto dst = Image::Make(odd ? h : w, odd ? w : h);

	float r, g, b, a;
	dst->Apply([&](int i, int j, float& r, float& g, float& b, float& a)
		{
			if (rots90 == 0) img->Get(i, j, r, g, b, a);
			if (rots90 == 1) img->Get(w - 1 - j, i, r, g, b, a);
			if (rots90 == 2) img->Get(w - 1 - i, h - 1 - j, r, g, b, a);
			if (rots90 == 3) img->Get(j, h - 1 - i, r, g, b, a);
		}
	);
	return dst;
}


struct Filter {};

ImagePtr xshear(const ImagePtr& img, float s, Filter kind/*='linear'*/)
{
	auto w = img->width;
	auto offset = (int)(ceil(s * img->width)); // todo - can fail?
	// this tweaks the placement so the center always ends up in the actual center
	auto shift = (offset - s * img->width) / 2.0f;
	auto newsize_w = img->width;
	auto newsize_h = abs(offset) + img->height;
	if (offset < 0)
		offset = abs(offset);
	else
		offset = 0;
	auto out = Image::Make(newsize_w, newsize_h);
	for (int row = 0; row < img->width; ++row)
	{
		throw image_exception("not impl");
		//		// xin,yin=nanfilter(np.asarray(range(img.shape[1])),img[row,:])
		//		auto f = interp1d(np.asarray(range(img.shape[1])), img[row, :], kind = kind, bounds_error = False, assume_sorted = True)
		//		auto pts = np.asarray(range(newsize[1])) - shift - offset - s * np.float(row)
		//		out[row, :] = map(clamp, f(pts))
	}
	return out;
}

ImagePtr yshear(const ImagePtr& img, float s, Filter kind/* = 'linear'*/)
{
	//inefficient, but avoids code duplication - todo - make fast
	auto tmp1 = RotateOrtho(img, 1);
	auto tmp2 = xshear(tmp1, s, kind);
	return RotateOrtho(tmp2, -1);
}

ImagePtr unser_rotate(const ImagePtr& img, float angle, Filter kind)
{
	auto xs = tanf(angle / 2.0f);
	auto ys = sinf(angle);
	auto tmp1 = xshear(img, xs, kind);
	auto tmp2 = yshear(tmp1, ys, kind);
	auto tmp3 = xshear(tmp2, xs, kind);
	auto [a, b] = tmp3->Size();
	auto [c, d] = img->Size();
	return Draw::Crop(tmp3,
		(a - c) / 2, (b - d) / 2,
		(a - c) / 2 + c, (b - d) / 2 + d
	);
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
ImagePtr ShearY(const ImagePtr& img, float shear)
{
	// todo- make faster
	auto rot1 = RotateOrtho(img, -1);
	auto sh = ShearX(rot1, shear);
	return RotateOrtho(sh, +1);
}

ImagePtr Rotate3Shear(const ImagePtr& img, float angleInRadians, bool expand)
{
	// jupyter nb at https://github.com/str4w/ImageRotation_PWL/blob/master/UnserRotation.ipynb

	// shear in x by -tan(angle/2) (explodes!)
	// shear in y by sin(angle)
	// shear in x by -tan(angle/2) (again! explodes!)

	// equivalent angle in [0,2pi)

	// also implement FULLY REVERSIBLE IMAGE ROTATION BY 1-D FILTERING - 2008 - better filter, works nicely

	auto a1 = fmodf(angleInRadians, pi * 2); // -2pi to +2pi
	while (a1 < 0) a1 += 2 * pi; // note 'if' may not be enough? (todo - check carefully)
	assert(a1 >= 0);
	auto q1 = fmodf(a1, pi / 2); // excess into quadrant
	assert(q1 >= 0);
	int rot90 = roundf((a1 - q1) / (pi / 2)); // num of rots 90

	// final fixup: instead of 0-pi/2 rotation, make it [-pi/4,pi/4]
	if (q1 > pi / 4)
	{
		rot90++;
		q1 -= pi / 2;
	}
	const ImagePtr src = (rot90 > 0) ? RotateOrtho(img, rot90) : img;

	auto shearX = -tanf(q1 / 2);
	auto shearY = sinf(q1);
	auto s1 = ShearX(src, shearX);
	auto s2 = ShearY(s1, -shearY);
	auto s3 = ShearX(s2, shearX);
	if (expand)
		return s3;
	return Draw::CropToBounds(s3); // todo - image was too large, clip more intelligently
}


ImagePtr Rotate(const ImagePtr& img, float angleInRadians, RotMethod method, bool expand)
{
	switch (method)
	{
	case RotMethod::DestNearest:
		return RotateDest(img, angleInRadians, expand, InterpNN);
	case RotMethod::DestBilinear:
		return RotateDest(img, angleInRadians, expand, InterpBilin);
	case RotMethod::DestShiftedBilinearA:
		return RotateDest(img, angleInRadians, expand, InterpShiftedBilinA);
	case RotMethod::DestShiftedBilinearB:
		return RotateDest(img, angleInRadians, expand, InterpShiftedBilinB);
	case RotMethod::DestBicubic:
		return RotateDest(img, angleInRadians, expand, InterpBicub);
	case RotMethod::DestDCT:
		return RotateDest(img, angleInRadians, expand, InterpDCT);
	case RotMethod::ShearShift:
		return Rotate3Shear(img, angleInRadians, expand);
	default:
		throw image_exception("Unknown rotation type");
	}
}

#endif