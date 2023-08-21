#pragma once
#include <string>
#include "State.h"
#include "ImgLib/ImgLib.h"
#include "fmt/fmt/format.h"
using namespace std;

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
		const int w = wi, h = hi; // clang has issue with auto vars...

		if (w2 == 0)
		{
			// w2/h2 = w/h
			w2 = static_cast<int>(round((double)h2 * w / h));
		}
		else if (h2 == 0)
		{
			// w2/h2 = w/h
			h2 = static_cast<int>(round((double)w2 * h / w));
		}
	}
	else if (args == "resize%")
	{
		const auto pct = s.Pop<double>();
		img = s.Pop<ImagePtr>();
		auto [w, h] = img->Size();
		w2 = static_cast<int>(round(w * pct / 100.0));
		h2 = static_cast<int>(round(h * pct / 100.0));
	}
	else if (args == "resize*")
	{
		const auto m = s.Pop<double>();
		img = s.Pop<ImagePtr>();
		auto [w, h] = img->Size();
		w2 = static_cast<int>(round(w * m));
		h2 = static_cast<int>(round(h * m));
	}
	else
		throw runtime_error(fmt::format("Unknown resize arg {}", method));
	auto [w, h] = img->Size();
	if (s.verbosity >= 1)
		cout << fmt::format("Resizing {}, {}x{} -> {}x{}, ", method, w, h, w2, h2);
	const Timer timer;
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
	else if (method == "lanczos2r")
		img = ResizeLanczosRadial(img, w2, h2, 2.0);
	else if (method == "lanczos3r")
		img = ResizeLanczosRadial(img, w2, h2, 3.0);
	else if (method == "lanczos4r")
		img = ResizeLanczosRadial(img, w2, h2, 4.0);
	else
		throw runtime_error(fmt::format("Unknown resize method {}", method));
	const auto elapsed = timer.get_elapsed_time();
	if (s.verbosity >= 1)
		cout << Timer::format_us(elapsed) << endl;
	s.Push(img);
};

/* -------------------- Misc image ops -----------------------*/

void GaussianBlur(State& s, const string& args)
{ // todo - use separable, make classes to abstract out kernels as templates
	const auto radius = s.Pop<double>();
	const auto src = s.Pop<ImagePtr>();
	auto [w1, h1] = src->Size();
	const int w = w1, h = h1; // clang glitch
	const auto dst = Image::Make(w, h);

	const double sigma = 1.0; // todo - base on kernel length?
	const double eConst = 2 * sigma * sigma;
	const double gConst = 1.0 / (eConst * std::numbers::pi);

	dst->Apply([=](int i, int j) {
		Color color(0, 0, 0, 0);
		double sum = 0;
		for (int dj = -static_cast<int>(floor(radius)); dj <= static_cast<int>(ceil(radius)); ++dj)
			for (int di = -static_cast<int>(floor(radius)); di <= static_cast<int>(ceil(radius)); ++di)
			{
				const auto d2 = di * di + dj * dj;
				if (d2 >= radius) continue;
				const auto weight = gConst * exp(d2 * eConst);
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
	const auto a = s.PopInt();
	const auto b = s.PopInt();
	const auto g = s.PopInt();
	const auto r = s.PopInt();
	auto right = s.PopInt();
	auto left = s.PopInt();
	auto bottom = s.PopInt();
	auto top = s.PopInt();
	const auto img = s.Pop<ImagePtr>();
	Color c(
		clamp(r, 0, 255) / 255.0,
		clamp(g, 0, 255) / 255.0,
		clamp(b, 0, 255) / 255.0,
		clamp(a, 0, 255) / 255.0
	);

	right = max(right, 0);
	left = max(left, 0);
	top = max(top, 0);
	bottom = max(bottom, 0);

	auto [w, h] = img->Size();
	const auto dst = Image::Make(left + w + right, top + h + bottom);
	dst->Apply([&](int i, int j) {
		if (i < left || j < top)
		{

		}
		return img->Get(i - left, j - top);
		});
	s.Push(dst);
}

void FlipImage(State& s, const string& args)
{
	const auto img = s.Pop<ImagePtr>();
	auto [w1, h1] = img->Size();
	int w = w1, h = h1; // todo - needed for clang issues on using w1,h1 in lambdas
	const auto dst = Image::Make(w, h);
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
	const auto src = s.Pop<ImagePtr>();
	auto [w, h] = src->Size();
	int w2 = x2 - x1 + 1, h2 = y2 - y1 + 1;
	const auto dst = Image::Make(w2, h2);
	for (int j = 0; j < h2; ++j)
		for (int i = 0; i < w2; ++i)
			dst->Set(i, j, src->Get(i + x1, j + y1));
	s.Push(dst);
}


void ImageOp(State& s, const string& args)
{
	if (args == "read")
	{
		string filename = s.Pop<string>();
		if (s.verbosity >= 1)
			cout << "Reading " << filename << endl;
		auto image = Image::Make(filename);
		s.Push(image);
	}
	else if (args == "write")
	{
		vector<uint8_t> data;
		auto filename = s.Pop<string>();
		auto img = s.Pop<ImagePtr>();
		cout << "Writing " << filename << endl;
		img->Save(filename);
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
		auto img = Image::Make(w, h);
		img->Apply([=](int i, int j) { return color; });
		s.Push(img);
	}
	else if (args == "i->f")
	{
		//		{"i->f", "f1 f2 .. fn n -> i1 i2 .. in, converts n values in 0-1 to n values in 0-255, useful for colors", ImageOp},
		int n = s.PopInt();
		vector<int> v;
		for (int i = 0; i < n; ++i)
			v.push_back(s.PopInt());
		for (int i = n - 1; i >= 0; --i)
			s.Push(Image::iToF64(v[i]));
	}
	else if (args == "f->i")
	{
		//		{ "f->i","i1 i2 .. in n -> f1 f2 .. fn, converts n values in 0.255 to n values in 0-1, useful for colors",ImageOp },
		int n = s.PopInt();
		vector<double> v;
		for (int i = 0; i < n; ++i)
			v.push_back(s.Pop<double>());
		for (int i = n - 1; i >= 0; --i)
			s.Push((int)(Image::f64ToI(v[i])));
	}
	else if (args == "blit" || args == "blitc" || args == "blitr" || args == "blitover")
	{
		int dx = 0, dy = 0, x1 = 0, y1 = 0, w = 0, h = 0;
		bool alphaBlend = false;
		ImagePtr overImage{ nullptr }, underImage{ nullptr };
		if (args == "blit")
		{
			//{"blit", "src dst -> dst', copy pixels from src to dst", ImageOp},
			underImage = s.Pop<ImagePtr>();
			overImage = s.Pop<ImagePtr>();
			auto [w1, h1] = overImage->Size();
			w = w1; h = h1;
		}
		else if (args == "blitc")
		{
			//{ "blitc", "src dst dx dy -> dst' copy src pixels to dst, placing dest corner at dx dy", ImageOp },
			dy = s.PopInt();
			dx = s.PopInt();
			underImage = s.Pop<ImagePtr>();
			overImage = s.Pop<ImagePtr>();
			auto [w1, h1] = overImage->Size();
			w = w1; h = h1;
		}
		else if (args == "blitr")
		{
			//{ "blitr", "src x1 y1 w h dst dx dy -> dst', copy rect from src x1 y1 w h to dst at dx dy", ImageOp },
			dy = s.PopInt();
			dx = s.PopInt();
			underImage = s.Pop<ImagePtr>();

			h = s.PopInt();
			w = s.PopInt();
			y1 = s.PopInt();
			x1 = s.PopInt();
			overImage = s.Pop<ImagePtr>();
		}
		else if (args == "blitover")
		{
			//{ "blitover", "src dst -> dst dx dy', alpha blend src OVER dst, at dx dy", ImageOp },
			dy = s.PopInt();
			dx = s.PopInt();
			underImage = s.Pop<ImagePtr>();
			overImage = s.Pop<ImagePtr>();
			auto [w1, h1] = overImage->Size();
			w = w1; h = h1;
			alphaBlend = true;

		}

		auto result = Blit(underImage, dx, dy, overImage, x1, y1, w, h, alphaBlend);

		s.Push(result);
	}
	else if (args == "alpha*" || args == "alpha/")
	{
		auto src = s.Pop<ImagePtr>();
		// clone
		auto img = Image::Make(src); // clones
		if (args == "alpha*")
		{
			AlphaCorrect(img, true);
		}
		else 
		{
			AlphaCorrect(img, false);
		}
		s.Push(img);
	}
	else if (args == "getpixel")
	{
		const auto j = s.PopInt();
		const auto i = s.PopInt();
		const auto img = s.Pop<ImagePtr>();
		const auto c = img->Get(i, j);
		s.Push(img);
		s.Push(c.r);
		s.Push(c.g);
		s.Push(c.b);
		s.Push(c.a);
	}
	else if (args == "setpixel")
	{
		// { "setpixel", "img i j r g b a -> img, writes pixel", PixelOp },
		const auto a = s.Pop<double>();
		const auto b = s.Pop<double>();
		const auto g = s.Pop<double>();
		const auto r = s.Pop<double>();
		const auto j = s.PopInt();
		const auto i = s.PopInt();
		const auto img = s.Pop<ImagePtr>();
		const Color c(r, g, b, a);
		img->Set(i, j, c);
		s.Push(img);
	}

	else if (args == "boundary")
	{
		// {"boundary", "img mode -> img', set sample boundary mode to clamp, reflect, reverse, tile", ImageOp },
		auto mode = s.Pop<string>();
		BoundaryMode bmode;
		if (mode == "color")
		{
			bmode.mode = BoundaryMode::Mode::Color;
			double a = s.Pop<double>();
			double b = s.Pop<double>();
			double g = s.Pop<double>();
			double r = s.Pop<double>();
			bmode.color = Color(r, g, b, a);
		}
		else if (mode == "clamped")
			bmode.mode = BoundaryMode::Mode::Clamped;
		else if (mode == "reflect")
			bmode.mode = BoundaryMode::Mode::Reflect;
		else if (mode == "reverse")
			bmode.mode = BoundaryMode::Mode::Reverse;
		else if (mode == "tile")
			bmode.mode = BoundaryMode::Mode::Tile;
		else throw runtime_error("Unknown boundary mode");

		auto img = s.Pop<ImagePtr>();

		img->boundaryMode = bmode;

		s.Push(img);
	}
	else if (args == "maxc")
	{
		const auto img = s.Pop<ImagePtr>();
		double maxc = 0;
		img->Apply([&](int i, int j) {
			const auto c = img->Get(i, j);
			maxc = max(c.r, maxc);
			maxc = max(c.g, maxc);
			maxc = max(c.b, maxc);
			return c;
			});
		if (s.verbosity >= 1)
			cout << fmt::format("maxc {:0.3f}\n", maxc);
		s.Push(maxc);
	}
	else if (args == "error")
	{

		auto method = s.Pop<string>();
		const auto img2 = s.Pop<ImagePtr>();
		const auto img1 = s.Pop<ImagePtr>();
		auto [w1, h1] = img1->Size();
		auto [w2, h2] = img2->Size();
		if (w1 != w2 || h1 != h2)
			throw runtime_error(fmt::format("Error size mismatch {}x{} vs {}x{}", w1, h1, w2, h2));
		double err = 0;
		if (method == "mse")
			err = MetricMSE(img1, img2);
		else if (method == "psnr")
			err = MetricPSNR(img1, img2);
		else if (method == "ssim")
			err = MetricSSIM(img1, img2);
		else
			cout << fmt::format("Error - unknown metric {} {}", method, err) << endl;
		if (s.verbosity >= 1)
			cout << fmt::format("{}: {:0.3f}\n", method, err);
		s.Push(img1);
		s.Push(img2);
		s.Push(err);
	}
	else if (args == "colorspace")
	{
		auto method = s.Pop<string>();
		const auto img1 = s.Pop<ImagePtr>();
		auto [w, h] = img1->Size();
		const auto img2 = Image::Make(w, h); // don't overwrite original!
		if (method == "linear")
			img2->Apply([&](int i, int j) {auto c = img1->Get(i, j); c.ApplyRGB(ToLinear); return c; });
		else if (method == "sRGB")
			img2->Apply([&](int i, int j) {auto c = img1->Get(i, j); c.ApplyRGB(FromLinear); return c; });
		else if (method == "RGB")
			img2->Apply([&](int i, int j) {return RGB(img1->Get(i, j)); });
		else if (method == "YCbCr")
			img2->Apply([&](int i, int j) {return YCbCr(img1->Get(i, j)); });
		else
			throw runtime_error(fmt::format("Unknown color operation {}", method));
		s.Push(img2);
	}
	else if (args == "resize" || args == "resize%" || args == "resize*")
	{
		ResizeImage(s, args);
	}
	else throw runtime_error("Unknown image op");
}


void RotateImage(State& s, const string& args)
{
	// {"rotate", "img angle filter -> img', rotate image by angle degrees using filter (see resize)", RotateImage},
	const auto filter = s.Pop<string>();
	const auto degrees = s.Pop<double>();
	auto img = s.Pop<ImagePtr>();

	const auto radians = degrees * numbers::pi / 180.0;

	if (filter == "nn")
	{
		img = RotateDest(img, radians, false, InterpNN);
	}
	else if (filter == "bilinear")
	{
		img = RotateDest(img, radians, false, InterpBilinear);
	}
	else if (filter == "bicubic")
	{
		img = RotateDest(img, radians, false, InterpBicubic);
	}
	else
		throw runtime_error("Unsupported rotate filter");


	s.Push(img);
}

void ShiftImage(State& s, const string& args)
{
	const auto filter = s.Pop<string>();
	const auto dy = s.Pop<double>();
	const auto dx = s.Pop<double>();
	auto img = s.Pop<ImagePtr>();
	if (filter == "nn")
	{
		img = ShiftImage2(img, dx, dy, InterpNN);
	}
	else if (filter == "bilinear")
	{
		img = ShiftImage2(img, dx, dy, InterpBilinear);

	}
	else if (filter == "bicubic")
	{
		img = ShiftImage2(img, dx, dy, InterpBicubic);
	}
	else throw runtime_error("Unsupported filter in ShiftImage");

	s.Push(img);

}


void DrawOp(State& s, const string& args)
{
	if (args == "line")
	{
		//	{"line", "img x1 y1 x2 y2 r g b a -> img with line", DrawOp},
		const auto a = s.Pop<double>();
		const auto b = s.Pop<double>();
		const auto g = s.Pop<double>();
		const auto r = s.Pop<double>();
		const auto y2 = s.PopInt();
		const auto x2 = s.PopInt();
		const auto y1 = s.PopInt();
		const auto x1 = s.PopInt();
		const auto img = s.Pop<ImagePtr>();
		DrawLine(img, x1, y1, x2, y2, Color(r, g, b, a));
		s.Push(img);
	}
	else if (args == "circle" || args == "circlef")
	{
		//	{ "circle",  "img x1 y1 radius r g b a -> img with circle",DrawOp },
		//	{ "circlef", "img x1 y1 radius r g b a -> img with filled circle",DrawOp },
		const auto a = s.Pop<double>();
		const auto b = s.Pop<double>();
		const auto g = s.Pop<double>();
		const auto r = s.Pop<double>();
		const auto radius = s.PopInt();
		const auto y1 = s.PopInt();
		const auto x1 = s.PopInt();
		const auto img = s.Pop<ImagePtr>();
		DrawCircle(img, x1, y1, radius, Color(r, g, b, a), args == "circlef");
		s.Push(img);

	}
	else if (args == "rect" || args == "rectf")
	{
		//	{ "rect",    "img x1 y1 x2 y2 r g b a -> img with rectangle",DrawOp },
		//	{ "rectf",   "img x1 y1 x2 y2 r g b a -> img with filled rectangle",DrawOp },
		const auto a = s.Pop<double>();
		const auto b = s.Pop<double>();
		const auto g = s.Pop<double>();
		const auto r = s.Pop<double>();
		const auto y2 = s.PopInt();
		const auto x2 = s.PopInt();
		const auto y1 = s.PopInt();
		const auto x1 = s.PopInt();
		const auto img = s.Pop<ImagePtr>();
		DrawRect(img, x1, y1, x2, y2, Color(r, g, b, a), args == "rectf");
		s.Push(img);
	}
	else if (args == "text")
	{
		//	{ "text",    "img x1 y1 r g b a text 0 m -> img x2 y2, draws text in font (always 0), pixel size m, returns img and final position",DrawOp },
		const auto m = s.PopInt();
		auto font = s.PopInt();
		const auto text = s.Pop<string>();
		const auto a = s.Pop<double>();
		const auto b = s.Pop<double>();
		const auto g = s.Pop<double>();
		const auto r = s.Pop<double>();
		const auto y1 = s.PopInt();
		const auto x1 = s.PopInt();
		const auto img = s.Pop<ImagePtr>();
		int dx, dy;
		DrawText(img, x1, y1, Color(r, g, b, a), text, m, dx, dy);
		s.Push(img);
		s.Push(dx);
		s.Push(dy);
	}
	else if (args == "gaussian")
	{
		GaussianBlur(s, args);
	}
	else if (args == "rotate") {
		RotateImage(s, args);
	}
	else if (args == "shift") {
		ShiftImage(s, args);
	}
	else if (args == "crop") {
		CropImage(s, args);
	}
	else if (args == "pad") {
		PadImage(s, args);
	}
	else if (args == "flipx" || args == "flipy") {
		FlipImage(s, args);
	}
	else throw runtime_error("Unknown draw op");

}
