#pragma once
#include <algorithm>
#include <functional>
#include <stdexcept>

/*---------------------- Color -----------------------------*/
// color treated exactly as a 4 vector, which works nicely with pre-multiplied alpha
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
		throw std::runtime_error("Invalid color channel index");
	}
	double& operator[](int index)
	{
		if (index == 0) return r;
		if (index == 1) return g;
		if (index == 2) return b;
		if (index == 3) return a;
		throw std::runtime_error("Invalid color channel index");
	}
	void ApplyRGB(const std::function<double(double)>& f)
	{
		r = f(r);
		g = f(g);
		b = f(b);
	}
	void ApplyRGBA(const std::function<double(double)>& f)
	{
		r = f(r);
		g = f(g);
		b = f(b);
		a = f(a);
	}
	Color& operator +=(const Color& c)
	{
		r += c.r;
		g += c.g;
		b += c.b;
		a += c.a;
		return *this;
	}
	Color& operator /=(double v)
	{
		r /= v;
		g /= v;
		b /= v;
		a /= v;
		return *this;
	}
	// clamp all internals to 0,1
	void Clamp()
	{
		r = std::clamp(r, 0.0, 1.0);
		g = std::clamp(g, 0.0, 1.0);
		b = std::clamp(b, 0.0, 1.0);
		a = std::clamp(a, 0.0, 1.0);
	}

	// return a new color that is this one with premultiplied alpha
	Color ToPremultipliedAlpha() const
	{
		if (std::abs(a) < 1e-20) return Color(0,0,0,0); // zero it all out
		return Color(r * a, g * a, b * a, a);
	}
	Color FromPremultipliedAlpha() const
	{
		if (std::abs(a) < 1e-20) return Color(0, 0, 0, 0); // zero it all out
		return Color(r / a, g / a, b / a, a);
	}

};

inline Color operator *(double w, const Color& c)
{
	return Color(w * c.r, w * c.g, w * c.b, w * c.a);
}

inline Color operator *(const Color& c, double w)
{
	return w * c;
}

inline Color operator +(const Color& c1, const Color& c2)
{
	return Color(c1.r + c2.r, c1.g + c2.g, c1.b + c2.b, c1.a + c2.a);
}