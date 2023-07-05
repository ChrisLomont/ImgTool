#pragma once
#include <algorithm>
#include <string>
#include "State.h"

using namespace std; // todo - remove

double ToLinear(double v)
{
	if (v <= 0.04045) return clamp(v / 12.92, 0.0, 1.0);
	return clamp(pow((v + 0.055) / 1.055, 2.4), 0.0, 1.0);
}
double FromLinear(double v)
{
	if (v <= 0.0031308) return clamp(12.92 * v, 0.0, 1.0);
	return clamp(1.055 * pow(v, 1.0/2.4) - 0.055, 0.0, 1.0);
}

// these done whether linear or not
// use must know what to do!
Color YCbCr(const Color& rgb)
{ // https://www.mir.com/DMG/ycbcr.html
	const auto delta = 0.5; // recenter to 0-1
	const auto& c = rgb;

	auto y = 0.299 * c.r + 0.587 * c.g + 0.114 * c.b;
	auto cb = -0.168736 * c.r + -0.331264 * c.g + 0.500 * c.b + delta;
	auto cr = 0.500 * c.r + -0.418688 * c.g + -0.081312 * c.b + delta;
	return Color(y, cb, cr, c.a);
}
Color RGB(const Color & yCbCr)
{
	const auto delta = 0.5; // recenter to 0-1

	const auto y = yCbCr.r;
	const auto cb = yCbCr.g - delta;
	const auto cr = yCbCr.b - delta;

	auto r = 1.0 * y + 0 * cb + 1.402 * cr;
	auto g = 1.0 * y + -0.344136 * cb + -0.714136 * cr;
	auto b = 1.0 * y + 1.772 * cb + 0 * cr;
	return Color(r, g, b, yCbCr.a);
}

void ColorTransform(State& s, const string& args)
{
	auto method = s.Pop<string>();
	auto img1 = s.Pop<ImagePtr>();
	auto [w, h] = img1->Size(); 
	auto img2 = make_shared<Image>(w,h); // don't overwrite original!
	if (method == "linear")
		img2->Apply([&](int i, int j) {auto c = img1->Get(i, j); c.ApplyRGB(ToLinear); return c; });
	else if (method == "sRGB")
		img2->Apply([&](int i, int j) {auto c = img1->Get(i, j); c.ApplyRGB(FromLinear); return c; });
	else if (method == "RGB")
		img2->Apply([&](int i, int j) {return RGB(img1->Get(i, j)); });
	else if (method == "YCbCr")
		img2->Apply([&](int i, int j) {return YCbCr(img1->Get(i, j)); });
	else
		throw runtime_error(format("Unknown color operation {}", method));
	s.Push(img2);
}