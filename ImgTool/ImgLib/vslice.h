#pragma once

#include <vector>

// non-owning view of stride sized items
// todo - implement out of bounds: clamp, default val, reflection, etc. this should be abstract, trait, etc.
// todo - 
template<typename T>
struct vslice
{
	T* base;
	size_t stride, offset, num;

	// make from vector - vector must outlive slice
	vslice(std::vector<T>& v)
		: base(v.data())
		, stride(1)
		, offset(0)
		, num(v.size())
	{
	}

	vslice(
		T* base,
		std::size_t offset,
		std::size_t stride,
		std::size_t num
	) : base(base), stride(stride), offset(offset), num(num)
	{
	}

	T& operator[](std::size_t index)
	{
		Check(index);
		return *(base + offset + index * stride);
	}
	T const& operator[](std::size_t index) const
	{
		Check(index);
		return *(base + offset + index * stride);
	}

	void Check(std::size_t index) const
	{
		if (index < 0 || num <= index)
			throw std::exception("out of bounds vslice");
	}

};