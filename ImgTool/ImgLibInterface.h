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
	else if (method == "lanczos2r")
		img = ResizeLanczosRadial(img, w2, h2, 2.0);
	else if (method == "lanczos3r")
		img = ResizeLanczosRadial(img, w2, h2, 3.0);
	else if (method == "lanczos4r")
		img = ResizeLanczosRadial(img, w2, h2, 4.0);
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
	auto dst = Image::Make(w, h);

	const double sigma = 1.0; // todo - base on kernel length?
	const double eConst = 2 * sigma * sigma;
	const double gConst = 1.0 / (eConst * std::numbers::pi);

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

	right = max(right, 0);
	left = max(left, 0);
	top = max(top, 0);
	bottom = max(bottom, 0);

	auto [w, h] = img->Size();
	auto dst = make_shared<Image>(left + w + right, top + h + bottom);
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
	auto img = s.Pop<ImagePtr>();
	auto [w1, h1] = img->Size();
	int w = w1, h = h1; // todo - needed for clang issues on using w1,h1 in lambdas
	auto dst = make_shared<Image>(w, h);
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
			auto c = img->Get(i, j);
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
		auto img = make_shared<Image>(w, h);
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
	else if (args == "blit" || args == "blitc" || args == "blitr")
	{
		int dx = 0, dy = 0, x1 = 0, y1 = 0, w = 0, h = 0;
		ImagePtr src{ nullptr }, dst{ nullptr };
		if (args == "blit")
		{
			//{"blit", "dst src -> dst', copy pixels from src to dst", ImageOp},
			src = s.Pop<ImagePtr>();
			dst = s.Pop<ImagePtr>();
			auto [w1, h1] = src->Size();
			w = w1; h = h1;
		}
		else if (args == "blitc")
		{
			//{ "blitc", "dst dx dy src -> dst' copy src pixels to dst, placing dest corner at dx dy", ImageOp },
			src = s.Pop<ImagePtr>();
			dy = s.PopInt();
			dx = s.PopInt();
			dst = s.Pop<ImagePtr>();
			auto [w1, h1] = src->Size();
			w = w1; h = h1;
		}
		else if (args == "blitr")
		{
			//{ "blitr", "dst dx dy src x1 y1 w h  -> dst', copy rect from src x1 y1 w h to dst at dx dy", ImageOp },
			h = s.PopInt();
			w = s.PopInt();
			y1 = s.PopInt();
			x1 = s.PopInt();

			src = s.Pop<ImagePtr>();
			dy = s.PopInt();
			dx = s.PopInt();
			dst = s.Pop<ImagePtr>();
		}

		Blit(dst, dx, dy, src, x1, y1, w, h);

		s.Push(dst);
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

	else throw runtime_error("Unknown image op");
}


void RotateImage(State& s, const string& args)
{
	// {"rotate", "img angle filter -> img', rotate image by angle degrees using filter (see resize)", RotateImage},
	auto filter = s.Pop<string>();
	auto degrees = s.Pop<double>();
	auto img = s.Pop<ImagePtr>();

	auto radians = degrees * numbers::pi / 180.0;

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
	auto filter = s.Pop<string>();
	auto dy = s.Pop<double>();
	auto dx = s.Pop<double>();
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



void ColorTransform(State& s, const string& args)
{
	auto method = s.Pop<string>();
	auto img1 = s.Pop<ImagePtr>();
	auto [w, h] = img1->Size();
	auto img2 = make_shared<Image>(w, h); // don't overwrite original!
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



void DrawOp(State& s, const string& args)
{
	if (args == "line")
	{
		//	{"line", "img x1 y1 x2 y2 r g b a -> img with line", DrawOp},
		auto a = s.Pop<double>();
		auto b = s.Pop<double>();
		auto g = s.Pop<double>();
		auto r = s.Pop<double>();
		auto y2 = s.PopInt();
		auto x2 = s.PopInt();
		auto y1 = s.PopInt();
		auto x1 = s.PopInt();
		auto img = s.Pop<ImagePtr>();
		DrawLine(img, x1, y1, x2, y2, Color(r, g, b, a));
		s.Push(img);
	}
	else if (args == "circle" || args == "circlef")
	{
		//	{ "circle",  "img x1 y1 radius r g b a -> img with circle",DrawOp },
		//	{ "circlef", "img x1 y1 radius r g b a -> img with filled circle",DrawOp },
		auto a = s.Pop<double>();
		auto b = s.Pop<double>();
		auto g = s.Pop<double>();
		auto r = s.Pop<double>();
		auto radius = s.PopInt();
		auto y1 = s.PopInt();
		auto x1 = s.PopInt();
		auto img = s.Pop<ImagePtr>();
		DrawCircle(img, x1, y1, radius, Color(r, g, b, a), args == "circlef");
		s.Push(img);

	}
	else if (args == "rect" || args == "rectf")
	{
		//	{ "rect",    "img x1 y1 x2 y2 r g b a -> img with rectangle",DrawOp },
		//	{ "rectf",   "img x1 y1 x2 y2 r g b a -> img with filled rectangle",DrawOp },
		auto a = s.Pop<double>();
		auto b = s.Pop<double>();
		auto g = s.Pop<double>();
		auto r = s.Pop<double>();
		auto y2 = s.PopInt();
		auto x2 = s.PopInt();
		auto y1 = s.PopInt();
		auto x1 = s.PopInt();
		auto img = s.Pop<ImagePtr>();
		DrawRect(img, x1, y1, x2, y2, Color(r, g, b, a), args == "rectf");
		s.Push(img);
	}
	else if (args == "text")
	{
		//	{ "text",    "img x1 y1 r g b a text 0 m -> img x2 y2, draws text in font (always 0), pixel size m, returns img and final position",DrawOp },
		auto m = s.PopInt();
		auto font = s.PopInt();
		auto text = s.Pop<string>();
		auto a = s.Pop<double>();
		auto b = s.Pop<double>();
		auto g = s.Pop<double>();
		auto r = s.Pop<double>();
		auto y1 = s.PopInt();
		auto x1 = s.PopInt();
		auto img = s.Pop<ImagePtr>();
		int dx, dy;
		DrawText(img, x1, y1, Color(r, g, b, a), text, m, dx, dy);
		s.Push(img);
		s.Push(dx);
		s.Push(dy);
	}
	else throw runtime_error("Unknown draw op");

}
