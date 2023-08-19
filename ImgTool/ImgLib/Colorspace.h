#pragma once
#include <algorithm>
#include <memory>

#include "Color.h"

using namespace std; // todo - remove

struct ColorspaceDef;
using ColorspaceDefPtr = std::shared_ptr<ColorspaceDef>;
struct ColorspaceDef
{
//	todo; // name
	static ColorspaceDefPtr Get(const std::string & name)
	{
		todo;
	}
};


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

	const auto y = 0.299 * c.r + 0.587 * c.g + 0.114 * c.b;
	const auto cb = -0.168736 * c.r + -0.331264 * c.g + 0.500 * c.b + delta;
	const auto cr = 0.500 * c.r + -0.418688 * c.g + -0.081312 * c.b + delta;
	return Color(y, cb, cr, c.a);
}
Color RGB(const Color & yCbCr)
{
	const auto delta = 0.5; // recenter to 0-1

	const auto y = yCbCr.r;
	const auto cb = yCbCr.g - delta;
	const auto cr = yCbCr.b - delta;

	const auto r = 1.0 * y + 0 * cb + 1.402 * cr;
	const auto g = 1.0 * y + -0.344136 * cb + -0.714136 * cr;
	const auto b = 1.0 * y + 1.772 * cb + 0 * cr;
	return Color(r, g, b, yCbCr.a);
}