#pragma once
#include <vector>
#include <iostream>
#include <algorithm>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#ifdef _MSC_VER  // needed for visual c++
#define __STDC_LIB_EXT1__ // sprintf -> sprintf_s
#endif
#include "stb_image_write.h"

#include "Timer.h"
#include "State.h"
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
			constexpr float pi2 = numbers::pi * numbers::pi;
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

ImagePtr ResizeLanczosRadial(ImagePtr src, int w, int h, double a)
{
	throw runtime_error("radial lanczos not implemented");
	// todo
}

double LanczosFunc(double v, int a)
{
	v = abs(v);
	if (v >= a) return 0;
	if (v == 0.0) return 1;
	double px = numbers::pi * v;
	// todo - for near 0, use taylor series, compute a few and put here
	return a * sin(px) * sin(px / a) / (px * px);
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


void ResizeImage(State& s, const string& args)
{ 
	auto method = ToLower(s.Pop<string>());
	ImagePtr img{ nullptr };
	int w2 = 0, h2 = 0;
	if (args == "resize")
	{
		// compute output width, height
		h2 = s.PopInt();
		w2 = s.PopInt();

		img = s.Pop<ImagePtr>();
		auto [wi, hi] = img->Size();
		int w = wi, h = hi; // clang has issue with auto vars...

		if (w2 == 0)
		{
			// w2/h2 = w/h
			w2 = (int)round((double)h2 * w / h);
		}
		else if (h2 == 0)
		{
			// w2/h2 = w/h
			h2 = (int)round((double)w2 * h / w);
		}
	}
	else if (args == "resize%")
	{
		auto pct = s.Pop<double>();
		img = s.Pop<ImagePtr>();
		auto [w, h] = img->Size();
		w2 = (int)(round(w * pct / 100.0));
		h2 = (int)(round(h * pct / 100.0));
	}
	else if (args == "resize*")
	{
		auto m = s.Pop<double>();
		img = s.Pop<ImagePtr>();
		auto [w, h] = img->Size();
		w2 = (int)(round(w * m));
		h2 = (int)(round(h * m));
	}
	else
		throw runtime_error(fmt::format("Unknown resize arg {}", method));
	auto [w, h] = img->Size();
	if (s.verbosity >= 1)
		cout << fmt::format("Resizing {}, {}x{} -> {}x{}, ", method, w, h, w2, h2);
	Timer timer;
	if (method == "nn")
		img = ResizeNN(img, w2, h2);
	else if (method == "bilinear")
		img = ResizeBilinear(img, w2, h2);
	else if (method == "bicubic")
		img = ResizeBicubic(img, w2, h2);
	else if (method == "lanczos2")
	{
		img = ApplyFilter<Lanczos2>(img, w2, h2);
		img = ResizeLanczos(img, w2, h2, 2);
	}
	else if (method == "lanczos3")
	{
		//img = ApplyFilter<Lanczos3>(img, w2, h2);
		img = ResizeLanczos(img, w2, h2, 3);
	}
	else if (method == "lanczos4")
	{
		//img = ApplyFilter<Lanczos4>(img, w2, h2);
		img = ResizeLanczos(img, w2, h2, 4);
	}
	else if (method == "lanczosr3")
		img = ResizeLanczosRadial(img, w2, h2, 3.0);
	else
		throw runtime_error(fmt::format("Unknown resize method {}", method));
	auto elapsed = timer.get_elapsed_time();
	if (s.verbosity >= 1)
		cout << Timer::format_us(elapsed) << endl;
	s.Push(img);
};

/* -------------------- Misc image ops -----------------------*/

void GaussianBlur(State& s, const string& args)
{ // todo - use separable, make classes to abstract out kernels as templates
	auto radius = s.Pop<double>();
	auto src = s.Pop<ImagePtr>();
	auto [w1, h1] = src->Size();
	int w = w1, h = h1; // clang glitch
	auto dst = Image::Make(w,h);

	const double sigma = 1.0; // todo - base on kernel length?
	const double eConst = 2 * sigma * sigma;
	const double gConst = 1.0 / (eConst*  std::numbers::pi);
	const double eConst = 2 * sigma * sigma;

	dst->Apply([=](int i, int j) {
		Color color(0, 0, 0, 0);
		double sum = 0;
		for (int dj = -(int)floor(radius); dj <= (int)ceil(radius); ++dj)
			for (int di = -(int)floor(radius); di <= (int)ceil(radius); ++di)
			{
				auto d2 = di * di + dj * dj;
				if (d2 >= radius) continue;
				auto weight = gConst * exp(d2 * eConst);
				sum += weight;
				color += weight * dst->Get(i + di, j + dj);
			}
		color /= sum;
		color.a = 1.0; // todo - more principled?
		return color;
		}
	);
	s.Push(dst);

}

void PadImage(State& s, const string& args)
{
	auto a = s.PopInt();
	auto b = s.PopInt();
	auto g = s.PopInt();
	auto r = s.PopInt();
	auto right = s.PopInt();
	auto left = s.PopInt();
	auto bottom = s.PopInt();
	auto top = s.PopInt();
	auto img = s.Pop<ImagePtr>();
	Color c(
		clamp(r, 0, 255) / 255.0,
		clamp(g, 0, 255) / 255.0,
		clamp(b, 0, 255) / 255.0,
		clamp(a, 0, 255) / 255.0
	);

	right  = max(right,0);
	left   = max(right, 0);
	top    = max(right, 0);
	bottom = max(right, 0);
	
	auto [w, h] = img->Size();
	auto dst = make_shared<Image>(left + w + right, top + h + bottom);
	dst->Apply([&](int i, int j) {
		if (i < left || j < top)
		{

		}
		return img->Get(i-left,j-top);
		});
	s.Push(dst);
}

void FlipImage(State& s, const string& args)
{
	auto img = s.Pop<ImagePtr>();
	auto [w1, h1] = img->Size();
	int w = w1, h = h1; // todo - needed for clang issues on using w1,h1 in lambdas
	auto dst = make_shared<Image>(w,h);
	if (args == "flipx")
	{
		// todo - does double work, clean
		dst->Apply([&](int i, int j) {
			return img->Get(w - 1 - i, j);
			}
		);
	}
	else if (args == "flipy")
	{
		// todo - does double work, clean
		dst->Apply([&](int i, int j) {
			return img->Get(i, h - 1 - j);
			}
		);

	}
	else
		throw runtime_error("Invalid flip");
	s.Push(dst);
}

void CropImage(State& s, const string& args)
{
	auto y2 = s.PopInt();
	auto x2 = s.PopInt();
	auto y1 = s.PopInt();
	auto x1 = s.PopInt();
	if (x2 < x1) swap(x1, x2);
	if (y2 < y1) swap(y1, y2);
	auto src = s.Pop<ImagePtr>();
	auto [w, h] = src->Size();
	int w2 = x2 - x1 + 1, h2 = y2 - y1 + 1;
	auto dst = make_shared<Image>(w2, h2);
	for (int j = 0; j < h2; ++j)
		for (int i = 0; i < w2; ++i)
			dst->Set(i, j, src->Get(i + x1, j + y1));
	s.Push(dst);
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

// should convert to gray first
// be careful about gamma or linear!
void ImageError(State& s, const string& args)
{
	string msg{ "" };
	if (args == "maxc")
	{
		auto img = s.Pop<ImagePtr>();
		double maxc = 0;
		img->Apply([&](int i, int j) {
			auto c = img->Get(i,j);
			maxc = max(c.r, maxc);
			maxc = max(c.g, maxc);
			maxc = max(c.b, maxc);
			return c;
			});
		msg = fmt::format("maxc {:0.3f}\n", maxc);
		s.Push(maxc);
	}
	else
	{

		auto method = s.Pop<string>();
		auto img2 = s.Pop<ImagePtr>();
		auto img1 = s.Pop<ImagePtr>();
		auto [w1, h1] = img1->Size();
		auto [w2, h2] = img2->Size();
		if (w1 != w2 || h1 != h2)
			throw runtime_error(fmt::format("Error size mismatch {}x{} vs {}x{}",w1,h1,w2,h2));
		double err = 0;
		if (method == "mse")
			err = MetricMSE(img1, img2);
		else if (method == "psnr")
			err = MetricPSNR(img1, img2);
		else if (method == "ssim")
			err = MetricSSIM(img1, img2);
		else
			cout << fmt::format("Error - unknown metric {} {}", method, err) << endl;
		msg = fmt::format("{}: {:0.3f}\n", method, err);
		s.Push(img1);
		s.Push(img2);
		s.Push(err);
	}
	if (s.verbosity >= 1)
		cout << msg;
}
void ImageOp(State& s, const string& args)
{
	if (args == "read")
	{
		int w, h, n;
		int channels = 4; // request RGBA
		string filename = s.Pop<string>();
		if (s.verbosity >= 1)
			cout << "Reading " << filename << endl;
		unsigned char* data = stbi_load(filename.c_str(), &w, &h, &n, channels);
		if (data == nullptr)
		{
			throw runtime_error(fmt::format("Cannot load file {}\n", filename));
		}
		//cout << "Data " << data << endl;
		auto img = make_shared<Image>(w, h, channels, data);
		stbi_image_free(data);
		s.Push(img);
	}
	else if (args == "write")
	{
		vector<uint8_t> data;
		auto filename = s.Pop<string>();
		auto img = s.Pop<ImagePtr>();
		img->GetData(data);
		auto [w, h] = img->Size();
		cout << "Writing " << filename << endl;

		int error = stbi_write_png(filename.c_str(),
			w, h, img->channels, data.data(), w * img->channels);
		s.Push(img);
	}
	else if (args == "size")
	{

		auto img = s.Peek<ImagePtr>();
		auto [w, h] = img->Size();
		s.Push(w);
		s.Push(h);
	}
	else if (args == "image")
	{
		// 	{"image","w h r g b a -> image, makes image size w x h, color rgba in 0-1",ImageOp},
		double a = s.Pop<double>();
		double b = s.Pop<double>();
		double g = s.Pop<double>();
		double r = s.Pop<double>();
		int h = s.PopInt();
		int w = s.PopInt();
		Color color(r, g, b, a);
		auto img = make_shared<Image>(w,h);
		img->Apply([=](int i, int j) { return color; });
		s.Push(img);
	}
	else throw runtime_error("Unknown image op");
}



void RotateImage(State& s, const string& args)
{
	throw runtime_error("Rotate functions not implemented");
}

#if 0
// DDA - lomont derivation
 
  /// <summary>
		/// Iterate over each x,y pair on the line between (x1,y1) and (x2,y2)
		/// Does so symmetrically: either direction makes the same line
		/// pixel midpoints are rounded towards positive
		/// </summary>
		/// <param name="x1"></param>
		/// <param name="y1"></param>
		/// <param name="x2"></param>
		/// <param name="y2"></param>
		public static IEnumerable<(int, int)> Dim2(int x1, int y1, int x2, int y2)
		{
			/* Derivation: Chris Lomont
			 Do 4 cases to handle slopes -1 to 1, varied on dx, dy, on x major
			 Two need "<" error comparisons, two need "<=" for symmetry
			 For y major, swap meaning of x,y except on output
			 Rewrite terms to get cases to have same internal loops
			 Reverse sign of e to get all comparisons in same direction
			 offset <= case by one to get to < case (or vice versa)
			 Merge cases
			 */

			 // todo - clean and simplify more

var(dx, dy) = (x2 - x1, y2 - y1);

// swap x,y to make slope in [-1,1]
var swap = Abs(dy) > Abs(dx);
if (swap)
{
	(dx, dy) = (dy, dx);
	(x1, y1) = (y1, x1);
	(x2, y2) = (y2, x2);
}

var(sx, sy) = (Sign(dx), Sign(dy));
var n = Max(Abs(dx), Abs(dy)) + 1;

var(x, y) = (x1, y1);

// in x major, stores 2*dx* error of (y-yi). is 0 for first pixel
// NOTE: while merging cases, meaning of e gets negated for some
// in y major, roles and x and y swapped
var e = 0;

// multiplying dx2 by sx*sy lets the code cases below match better
var(dx2, dy2) = (2 * dx * sx * sy, 2 * dy);

// error comparison for de <= e cases
var ec = dx * (sx);

if (
	(0 < dx && 0 <= dy && Abs(dy) <= Abs(dx)) || // slope 0 <= m <= 1
	(dx < 0 && 0 <= dy && Abs(dy) <= Abs(dx)) // slope -1 <= m <= 0
	)

{  // slope 0 <= m <= 1
}
else if (
	// slope 0 <= m <= 1
	(dx < 0 && dy <= 0 && Abs(dy) <= Abs(dx))
	||
	// slope -1 <= m <= 0
	(0 < dx && dy <= 0 && Abs(dy) <= Abs(dx))
	)
{
	// reverse e meaning
	(dx2, dy2) = (-dx2, -dy2);

	// shift ec to change "<" to "<="
	ec += 1;
}

if (swap)
{
	for (var i = 0; i < n; ++i)
	{
		yield return new(y, x);
		x += sx; // move point in dx direction
		e += dy2; // updated error from move
		if (ec <= e) // is error too big?
		{
			y += sy;  // move point up/down
			e -= dx2; // error gets adjusted
		}
	}
}
else
{
	for (var i = 0; i < n; ++i)
	{
		yield return new(x, y);
		x += sx; // move point in dx direction
		e += dy2; // updated error from move
		if (ec <= e) // is error too big?
		{
			y += sy; // move point up/down
			e -= dx2; // error gets adjusted
		}
	}
}
		}


#endif