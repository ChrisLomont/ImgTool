#pragma once
#include <vector>
#include <functional>
#include <algorithm>
#include <memory>
#include "Boundary.h"
#include "Color.h"
#include "vslice.h"
#include "ImageFormat.h"

using namespace std; // todo - remove this

/*---------------------- Image -----------------------------*/
class Image;
using ImagePtr = shared_ptr<Image>;

class Image {
	vector<Color> data_;
	int w, h;

	// read file helper - read into image
	static void ReadHelper(const std::string& filename, Image& img);


public:
	// todo - link lomont blog onn the reasoning
	// // todo - see Lomont article on quantization
	static double iToF64(int v) { return v / 255.0; }
	static int f64ToI(double v) { return static_cast<int>(clamp(floor(v * 256.0), 0.0, 255.0)); }

	// todo - implement and leverage
	ImageFormat imageFormat;

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

	Image(const Image& image)
	{
		Resize(image.w, image.h);
		data_ = image.data_; // copies
	}
	Image(const std::string& filename)
	{
		ReadHelper(filename, *this);
	}
	//Image(int width = 1, int height = 1)
	//{
	//	Resize(width, height);
	//}
	Image& operator=(const Image& img)
	{
		Resize(img.w, img.h);
		data_ = img.data_; // does copy
		return *this;
	}
	void Resize(int width, int height)
	{
		this->w = width;
		this->h = height;
		data_.resize(width * height * channels);
	}
	// make copy
	static ImagePtr Make(const ImagePtr& img)
	{
		auto [w, h] = img->Size();
		auto dst = Make(w, h);
		dst->data_ = img->data_; // does copy
		return dst;
	}

	static ImagePtr Make(const std::string& filename)
	{
		return std::make_shared<Image>(filename);
	}

	static ImagePtr Make(int w, int h)
	{
		return std::make_shared<Image>(w, h);
	}

	// read file
	static ImagePtr Read(const std::string& filename);
	void Save(const std::string& filename);


	[[nodiscard]] pair<int, int> Size() const { return pair(w, h); }

	// modify in place by pixel
	//void Apply(std::function<void(float& r, float& g, float& b, float& a)> func);
	// modify in place by pixel and index
	//void Apply(std::function<void(int i, int j, float& r, float& g, float& b, float& a)> func);

	void Apply(const function<Color (const Color& c)>& colorMap)
	{
		for (auto j = 0; j < h; ++j)
			for (auto i = 0; i < w; ++i)
			{
				auto c1 = Get(i, j);
				auto c2 = colorMap(c1);
				Set(i, j, c2);
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
		return {colors, col, w, h};
	}
	// get row, modifiable
	// todo - make these understand boundary conditions also
	vslice<Color> Row(int row)
	{
		const auto colors = (Color*)data_.data();
		return {colors, row * w, 1, w};
	}

};
