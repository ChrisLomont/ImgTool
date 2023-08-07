#pragma once
#include <vector>
#include <functional>
#include <algorithm>
#include <memory>
#include "Boundary.h"
#include "Color.h"
#include "vslice.h"

using namespace std; // todo - remove this

/*---------------------- Image -----------------------------*/

class Image {
	vector<Color> data_;
	int w, h;


public:
	// todo - link lomont blog onn the reasoning
	// // todo - see Lomont article on quantization
	static double iToF64(int v) { return v / 255.0; }
	static int f64ToI(double v) { return static_cast<int>(clamp(floor(v * 256.0), 0.0, 255.0)); }

	static shared_ptr<Image> Make(int w, int h)
	{
		return make_shared<Image>(w,h);
	}

	[[nodiscard]] bool Legal(int i, int j) const { return 0 <= i && 0 <= j && i < w && j < h; }

	int channels;
	Image(int w, int h, int channels = 4, const unsigned char* data = nullptr)
		: w(w)
		, h(h)
		, channels(channels)
	{
		data_.resize(w * h);
		if (data != nullptr)
		{
			for (auto j = 0; j < h; ++j)
				for (auto i = 0; i < w; ++i)
				{
					int index = (i + j * w) * channels;
					const double r = iToF64(data[index++]);
					const double g = iToF64(data[index++]);
					const double b = iToF64(data[index++]);
					const double a = iToF64(data[index++]);
					Set(i, j, Color(r, g, b, a));
				}
		}
	}

	[[nodiscard]] pair<int, int> Size() const { return pair(w, h); }

	void Apply(const function<void(Color& c)>& colorMap)
	{
		for (auto j = 0; j < h; ++j)
			for (auto i = 0; i < w; ++i)
			{
				auto c = Get(i, j);
				colorMap(c);
				Set(i, j, c);
			}
	}
	void Apply(const function<Color(int i, int j)>& colorFunc)
	{
		for (auto j = 0; j < h; ++j)
			for (auto i = 0; i < w; ++i)
			{
				Set(i, j, colorFunc(i, j));
			}
	}
	
	// todo boundary implementation
	BoundaryMode boundaryMode;
	void Set(int i, int j, const Color& c) { if (Legal(i, j))data_[i + j * w] = c; }

	[[nodiscard]] Color Get(int i, int j) const {
		if (Legal(i, j))
			return data_[i + j * w];
		if (boundaryMode.mode == BoundaryMode::Mode::Color)
			return boundaryMode.color;
		i = BoundaryClamp(boundaryMode, i, 0, w);
		j = BoundaryClamp(boundaryMode, j, 0, h);
		return data_[i + j * w];
	}


	void GetData(vector<uint8_t>& data) const
	{
		data.resize(w * h * channels);
		for (auto j = 0; j < h; ++j)
			for (auto i = 0; i < w; ++i)
			{
				const auto c = Get(i, j);

				int index = (i + j * w) * channels;

				data[index++] = f64ToI(c.r);
				data[index++] = f64ToI(c.g);
				data[index++] = f64ToI(c.b);
				data[index++] = f64ToI(c.a);

			}
	}

	// get column, modifiable
	vslice<Color> Col(int col)
	{
		const auto colors = (Color*)data_.data();
		return vslice<Color>(colors, col, w, h);
	}
	vslice<Color> Row(int row)
	{
		const auto colors = (Color*)data_.data();
		return vslice<Color>(colors, row * w, 1, w);
	}

};
using ImagePtr = shared_ptr<Image>;