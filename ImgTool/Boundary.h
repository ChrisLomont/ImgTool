#pragma once
#include <algorithm>

// boundary conditions
enum class BoundaryMode
{
	Clamp = 0,
	Reflect = 1, // reflect edge, edge pixel only occurs once
	Reverse = 2,  // reflect edge, edge pixel occurs twice
	Tile = 3
};

// compute positive part of v%m for any v,m
int PositiveMod(int v, int m)
{
	return ((v % m) + m) % m;
}

// perform clamping
// usually used on integer coords [0,w) for an image
// min inclusive, max exclusive
int BoundaryClamp(BoundaryMode mode, int val, int min, int max)
{
	if (mode == BoundaryMode::Clamp)
		return std::clamp(val,min,max-1);
	else if (mode == BoundaryMode::Tile)
	{ 
		return PositiveMod(val, max - min) + min;
	}
	else if (mode == BoundaryMode::Reverse)
	{ // repeat edge item
		const auto w = max - min;
		const auto m = 2 * w;
		auto v1 = PositiveMod(val, m);
		if (v1 >= w)
			v1 = 2 * w - v1 - 1;
		return v1;
	}
	else if (mode == BoundaryMode::Reflect)
	{
		// don't repeat edge item
		const auto w = max - min - 1;
		const auto m = 2 * w;
		auto v1 = PositiveMod(val, m);
		if (v1 >= w)
			v1 = 2 * w - v1;
		return v1;
	}
	return (min + max) / 2;
}
