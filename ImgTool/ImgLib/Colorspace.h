#pragma once
#include <algorithm>
#include <memory>

#include "Color.h"


struct ColorspaceDef;
using ColorspaceDefPtr = std::shared_ptr<ColorspaceDef>;
struct ColorspaceDef
{
//	todo; // name
	static ColorspaceDefPtr Get(const std::string & name)
	{
		//todo;
		return std::make_shared<ColorspaceDef>();
	}
};

// srgb to linear
double ToLinear(double v);
// linear to srgb
double FromLinear(double v);


// these done whether linear or not
// use must know what to do!
Color YCbCr(const Color& rgb);

Color RGB(const Color& yCbCr);