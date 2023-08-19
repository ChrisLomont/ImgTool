#pragma once

#include "Colorspace.h"

// a place to store image format info
enum class AlphaFormat
{
	NoAlpha,
	AssociatedAlpha,   // same as premultiplied alpha
	UnassociatedAlpha, // is not premultiplied alpha
};
struct ImageFormat
{
	// meaning of any alpha channel
	AlphaFormat alphaFormat{AlphaFormat::UnassociatedAlpha};

	// is the surface in linear or gamma color
	bool linearColor{ false };

	ColorspaceDefPtr colorspace = ColorspaceDef::Get("sRGB"); // default here
	
};