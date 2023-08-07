// Lomont old rotate code
# pragma once
#if 0
#include "../Image/Image.h"

// rotate by given number of 90 degree rotations
ImagePtr RotateOrtho(const ImagePtr& img, int rots90, bool expand = true);

enum class RotMethod {
	// walk dest pixels, back project, and interpolate:
	DestNearest,
	DestBilinear,
	DestShiftedBilinearA,
	DestShiftedBilinearB,
	DestBicubic,
	DestDCT,

	// 3 shear, linear separable filters
	ShearShift, // shift pixels, lossless
	//Shear??
};

// rotate about center
ImagePtr Rotate(const ImagePtr& img, float angleInRadians, RotMethod method, bool expand = true);
#endif