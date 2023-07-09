#pragma once
#include <vector>
#include <iostream>
#include <algorithm>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define __STDC_LIB_EXT1__ // sprintf -> sprintf_s
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
	throw runtime_error("radial lanczos not impelmented");
	// todo
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
		auto [w, h] = img->Size();

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
		img = ApplyFilter<Lanczos2>(img, w2, h2);
	else if (method == "lanczos3")
		img = ApplyFilter<Lanczos3>(img, w2, h2);
	else if (method == "lanczos4")
		img = ApplyFilter<Lanczos4>(img, w2, h2);
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
{
	throw runtime_error("Gaussian not implemented");
	// todo;
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
	auto [w, h] = img->Size();
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

void ReadImage(State& s, const string& args)
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
void WriteImage(State& s, const string& args)
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

void RotateImage(State& s, const string& args)
{
	throw runtime_error("Rotate functions not implemented");
}