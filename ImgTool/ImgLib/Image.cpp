#include "Image.h"
#include "image_exception.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#ifdef _MSC_VER  // needed for visual c++
#define __STDC_LIB_EXT1__ // sprintf -> sprintf_s
#endif
#include <cctype>

#include "stb_image_write.h"

void Image::ReadHelper(const std::string& filename, Image& img)
{
	int w, h, n;
	int channels = 4; // request RGBA
	unsigned char* data = stbi_load(filename.c_str(), &w, &h, &n, channels);
	if (data == nullptr)
	{
		throw runtime_error(string("Cannot load file ") + filename);
	}
	img = Image(w, h, channels, data);
	stbi_image_free(data);
}


void Image::Save(const std::string& filename)
{
	int comp = 4; // channels
	int stride_in_bytes = w * comp;
	vector<uint8_t> img; // RGBA
	img.resize(w * h * 4);

	for (int j = 0; j < h; ++j)
		for (int i = 0; i < w; ++i)
		{
			auto color = Get(i, j);

			int ind = (i + j * w) * channels;
			img[ind + 0] = f64ToI(color.r);
			img[ind + 1] = f64ToI(color.g);
			img[ind + 2] = f64ToI(color.b);
			img[ind + 3] = f64ToI(color.a);
		}

	std::string str = filename;
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	int err = 0;
	if (str.ends_with(".png"))
	{
		err = stbi_write_png(filename.c_str(), w, h, comp,
			img.data(),
			stride_in_bytes
		);
	}
	if (str.ends_with(".jpg") || str.ends_with(".jpeg"))
	{
		err = stbi_write_jpg(filename.c_str(), w, h, comp,
			img.data(),
			100 // quality max
		);
	}
	if (str.ends_with(".bmp"))
	{
		err = stbi_write_bmp(filename.c_str(), w, h, comp,
			img.data()
		);
	}
	if (str.ends_with(".tga"))
	{
		err = stbi_write_tga(filename.c_str(), w, h, comp,
			img.data()
		);
	}
	if (str.ends_with(".hdr"))
	{
		//todo - fixerr = stbi_write_hdr(filename.c_str(), w, h, comp, data_.data());
		err = 0;
	}
	if (err == 0) throw image_exception("unsupported file save format");
}
