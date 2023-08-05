#pragma once
// basic drawing
#include <string>
#include <functional>
#include <cmath>
#include "Font9x16.h"
#include "Image.h"

using namespace std;

int sign(int val)
{
	if (val < 0) return -1;
	if (val > 0) return +1;
	return 0;
}

// Iterate over each x,y pair on the line between (x1,y1) and (x2,y2)
// Does so symmetrically: either direction makes the same line
// pixel midpoints are rounded towards positive
void DDA2(int x1, int y1, int x2, int y2, function<void(int i,int j)> processCoord)
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

	auto dx = x2 - x1, dy = y2 - y1;

	// swap x,y to make slope in [-1,1]
	auto swap2 = abs(dy) > abs(dx);
	if (swap2)
	{
		swap(dx, dy);
		swap(x1, y1);
		swap(x2, y2);
	}

	auto sx = sign(dx), sy = sign(dy);
	auto n = max(abs(dx), abs(dy)) + 1;

	auto x = x1, y = y1;

	// in x major, stores 2*dx* error of (y-yi). is 0 for first pixel
	// NOTE: while merging cases, meaning of e gets negated for some
	// in y major, roles and x and y swapped
	auto e = 0;

	// multiplying dx2 by sx*sy lets the code cases below match better
	auto dx2 = 2 * dx * sx * sy, dy2 = 2 * dy;

	// error comparison for de <= e cases
	auto ec = dx * (sx);

	if (
		(0 < dx && 0 <= dy && abs(dy) <= abs(dx)) || // slope 0 <= m <= 1
		(dx < 0 && 0 <= dy && abs(dy) <= abs(dx)) // slope -1 <= m <= 0
		)

	{  // slope 0 <= m <= 1
	}
	else if (
		// slope 0 <= m <= 1
		(dx < 0 && dy <= 0 && abs(dy) <= abs(dx))
		||
		// slope -1 <= m <= 0
		(0 < dx && dy <= 0 && abs(dy) <= abs(dx))
		)
	{
		// reverse e meaning
		dx2 = -dx2;
		dy2 = -dy2;

		// shift ec to change "<" to "<="
		ec += 1;
	}

	if (swap2)
	{
		for (auto i = 0; i < n; ++i)
		{
			processCoord(y, x);
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
		for (auto i = 0; i < n; ++i)
		{
			processCoord(x, y);
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

void DrawLine(ImagePtr img, int x1, int y1, int x2, int y2, const Color & color)
{
	DDA2(x1, y1, x2, y2, 
		[&](int i, int j) {
			img->Set(i, j, color);
		}
		);
}

void DrawRect(ImagePtr img, int x1, int y1, int x2, int y2, const Color& color, bool filled)
{
	if (filled)
	{
		// y1 < y2
		if (y1 > y2)
		{
			swap(x1, x2);
			swap(y1, y2);
		}
		for (int y = y1; y <= y2; ++y)
		{
			DrawLine(img, x1, y, x2, y, color);
		}
	}
	else {

		DrawLine(img, x1, y1, x2, y1, color);
		DrawLine(img, x2, y1, x2, y2, color);
		DrawLine(img, x2, y2, x1, y2, color);
		DrawLine(img, x1, y2, x1, y1, color);
	}
}

void DrawGlyph(ImagePtr img, char c, int x1, int y1, int mult, int & x2, int & y2, const Color &color)
{
	// todo - can vastly speed this up
	// 9x16 font, 18 bytes each
	// each pair a1,a2 is pixels
	const uint8_t* p = Font9x16 + c * 18;

	for (int j = 0; j < 16; ++j)
		for (int i = 0; i < 9; ++i)
		{
			int index = 2 * i + (j / 8);
			uint8_t b = p[index];
			int pixel = (b >> (j & 7)) & 1;
			if (pixel == 0) continue;

			for (int sy = 0; sy < mult; ++sy)
				for (int sx = 0; sx < mult; ++sx)
				{
					img->Set(x1 + i * mult + sx, y1 + j * mult + sy, color);
				}
		}
	x2 = x1 + 9*mult;
	y2 = y1 + 16*mult;
}

void DrawText(ImagePtr img, int x0, int y0, const Color& color, const string& text, int multiplier, int& x2, int& y2)
{
	auto x = x0, y = y0;
	for (auto c : text)
	{
		DrawGlyph(img, c, x, y, multiplier, x2, y2,color);
		x = x2 + multiplier; // spacing
	}
}

// Function to put pixels
// at subsequence points
void Draw8(ImagePtr img, const Color & color, int xc, int yc, int x, int y, bool filled)
{
	if (filled)
	{
		DrawLine(img, xc + x, yc + y, xc - x, yc + y, color);
		DrawLine(img, xc + x, yc - y, xc - x, yc - y, color);
		DrawLine(img, xc + y, yc + x, xc - y, yc + x, color);
		DrawLine(img, xc + y, yc - x, xc - y, yc - x, color);
	}
	else
	{
		img->Set(xc + x, yc + y, color);
		img->Set(xc - x, yc + y, color);

		img->Set(xc + x, yc - y, color);
		img->Set(xc - x, yc - y, color);

		img->Set(xc + y, yc + x, color);
		img->Set(xc - y, yc + x, color);

		img->Set(xc + y, yc - x, color);
		img->Set(xc - y, yc - x, color);
	}
}


void DrawCircle(ImagePtr img, int xc, int yc, int radius, const Color& color, bool filled)
{ // bresenham algo
	auto t1 = radius / 16;
	auto x = radius;
	auto y = 0;
	while (x >= y)
	{
		Draw8(img, color, xc, yc, x, y, filled);
		y = y + 1;
		t1 = t1 + y;
		auto t2 = t1 - x;
		if (t2 >= 0)
		{
			t1 = t2;
			x = x - 1;
		}
	}
}

