#pragma once
#include <vector>
#include <functional>
#include <algorithm>
#include <memory>
#include <stdexcept>

using namespace std; // todo - remove this

/*------------------- Color and Image ----------------------*/
struct Color
{
	double r, g, b, a;
	Color(double r = 1, double g = 0, double b = 1, double a = 1)
		: r(r)
		, g(g)
		, b(b)
		, a(a)
	{
	}
	double operator[](int index) const
	{
		if (index == 0) return r;
		if (index == 1) return g;
		if (index == 2) return b;
		if (index == 3) return a;
		throw runtime_error("Invalid color channel index");
	}
	double& operator[](int index)
	{
		if (index == 0) return r;
		if (index == 1) return g;
		if (index == 2) return b;
		if (index == 3) return a;
		throw runtime_error("Invalid color channel index");
	}
	void ApplyRGB(const function<double(double)>& f)
	{
		r = f(r);
		g = f(g);
		b = f(b);
	}
	Color& operator +=(const Color& c)
	{
		r += c.r;
		g += c.g;
		b += c.b;
		// todo - alpha belnding?
		return *this;
	}
	Color& operator /=(double v)
	{
		r /= v;
		g /= v;
		b /= v;
		// todo - alpha?
		return *this;
	}
};
Color operator *(double w, const Color& c)
{
	return Color(w * c.r, w * c.g, w * c.g, c.a);
}

class Image {
	vector<Color> data_;
	int w, h;


public:
	static double iToF64(int v) { return v / 255.0; }
	static int f64ToI(double v) { return (int)(clamp(floor(v * 256.0), 0.0, 255.0)); }

	static shared_ptr<Image> Make(int w, int h)
	{
		return make_shared<Image>(w,h);
	}
	bool Legal(int i, int j) const { return 0 <= i && 0 <= j && i < w && j < h; }

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
					double r = iToF64(data[index++]);
					double g = iToF64(data[index++]);
					double b = iToF64(data[index++]);
					double a = iToF64(data[index++]);
					Set(i, j, Color(r, g, b, a));
				}
		}
	}
	pair<int, int> Size() const { return pair(w, h); }

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

	void Set(int i, int j, const Color& c) { if (Legal(i, j))data_[i + j * w] = c; }
	Color Get(int i, int j) const { if (Legal(i, j)) return data_[i + j * w]; return Color(1, 0, 1, 1); }

	void GetData(vector<uint8_t>& data) const
	{
		data.resize(w * h * channels);
		for (auto j = 0; j < h; ++j)
			for (auto i = 0; i < w; ++i)
			{

				auto c = Get(i, j);

				int index = (i + j * w) * channels;

				data[index++] = f64ToI(c.r);
				data[index++] = f64ToI(c.g);
				data[index++] = f64ToI(c.b);
				data[index++] = f64ToI(c.a);

			}
	}
};
using ImagePtr = shared_ptr<Image>;