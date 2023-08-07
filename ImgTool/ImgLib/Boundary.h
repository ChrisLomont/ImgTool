#pragma once
#include <algorithm>
#include "Color.h"

struct BoundaryMode {
	// boundary conditions
	enum class Mode
	{
		Color = 0,   // use color below
		Clamped = 1,
		Reflect = 2, // reflect edge, edge pixel only occurs once
		Reverse = 3, // reflect edge, edge pixel occurs twice
		Tile = 4     // 
	};
	Mode mode{Mode::Reflect};
	// color if mode is Color
	Color color{0,0,0,0};
};

// compute positive part of v%m for any v,m
inline int PositiveMod(int v, int m)
{
	return ((v % m) + m) % m;
}

// perform clamping
// usually used on integer coords [0,w) for an image
// min inclusive, max exclusive
// do not pass Color mode into here, check outside?
inline int BoundaryClamp(const BoundaryMode & mode, int val, int min, int max)
{
	if (mode.mode == BoundaryMode::Mode::Clamped)
		return std::clamp(val,min,max-1);
	else if (mode.mode == BoundaryMode::Mode::Tile)
	{ 
		return PositiveMod(val, max - min) + min;
	}
	else if (mode.mode == BoundaryMode::Mode::Reverse)
	{ // repeat edge item
		const auto w = max - min;
		const auto m = 2 * w;
		auto v1 = PositiveMod(val, m);
		if (v1 >= w)
			v1 = 2 * w - v1 - 1;
		return v1;
	}
	else if (mode.mode == BoundaryMode::Mode::Reflect)
	{
		// don't repeat edge item
		const auto w = max - min - 1;
		const auto m = 2 * w;
		auto v1 = PositiveMod(val, m);
		if (v1 >= w)
			v1 = 2 * w - v1;
		return v1;
	}
	throw runtime_error("Invalid boundary mode");
}
