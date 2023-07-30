#pragma once

/*---------------------- Color -----------------------------*/
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
	// clamp all internals to 0,1
	void Clamp()
	{
		r = clamp(r, 0.0, 1.0);
		g = clamp(g, 0.0, 1.0);
		b = clamp(b, 0.0, 1.0);
		a = clamp(a, 0.0, 1.0);

	}
};
Color operator *(double w, const Color& c)
{
	return Color(w * c.r, w * c.g, w * c.b, c.a);
}
Color operator *(const Color& c, double w)
{
	return w * c;
}
Color operator +(const Color& c1, const Color& c2)
{
	return Color(c1.r + c2.r, c1.g + c2.g, c1.b + c2.b, (c1.a + c2.a) / 2);
}