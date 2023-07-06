#pragma once
#include <variant>
#include <vector>
#include <string>
#include <iostream>
#include "Utils.h"

#include "Image.h"

using namespace std; // todfo - remove

using Item = variant<string, ImagePtr, double, int>;

string FormatItem(Item & item, bool prefixType)
{
	if (holds_alternative<ImagePtr>(item))
	{
		auto img = get<ImagePtr>(item);
		auto [w, h] = img->Size();
		return prefixType
			? fmt::format("Image {} x {}", w, h)
			: fmt::format("<Image>");
	}
	else if (holds_alternative<string>(item))
	{
		return prefixType
			? fmt::format("string {}", get<string>(item))
			: get<string>(item);
	}
	else if (holds_alternative<int>(item))
	{
		return prefixType
			? fmt::format("int {}", get<int>(item))
			: fmt::format("{}", get<int>(item));

	}
	else if (holds_alternative<double>(item))
	{
		return prefixType
			? fmt::format("double {}", get<double>(item))
			: fmt::format("{}", get<double>(item));
	}
	else
		throw runtime_error(fmt::format("Unknown print type"));
}


// computation state
class Stack {
	vector<Item> stack;
	int nextopen = 0; // next one available
	void Resize() {
		while (stack.size() <= nextopen)
			stack.push_back("<todo>");
	}
	tuple<Item,Item> Pop2()
	{
		auto v1 = Pop();
		auto v2 = Pop();
		return { v1, v2 };
	}
	tuple<Item, Item, Item> Pop3()
	{
		auto v1 = Pop();
		auto v2 = Pop();
		auto v3 = Pop();
		return { v1, v2, v3 };
	}
	void Push(initializer_list<Item> items)
	{
		for (auto& i : items)
			Push(i);
	}
protected:
	vector<Item> PopN(int n)
	{
		vector<Item> items;
		for (int i = 0; i < n; ++i)
			items.push_back(Pop());
		return items;
	}
	// push, REVERSE ORDER!
	void PushN(const vector<Item> & items)
	{
		if (items.size() == 0) return;
		for (int i = 0; i < items.size(); ++i)
			Push(items[items.size() - 1 - i]);
	}
public:
	int Pos() const { return nextopen; }

	void StackOp(const string& arg)
	{
		if (arg == "dup")
		{
			Push(Peek(0));
		}
		else if (arg == "dup2")
		{
			auto [a, b] = Pop2();
			Push({b,a,b,a});
		}
		else if (arg == "dupn")
		{
			auto n = ParseInt(Pop<string>());
			for (int i = 0; i < n; ++i)
				Push(stack[nextopen - n]);
		}
		else if (arg == "drop")
		{
			Pop();
		}
		else if (arg == "drop2")
		{
			Pop2();
		}
		else if (arg == "dropn")
		{
			auto n = ParseInt(Pop<string>());
			for (int i = 0; i < n; ++i)
				Pop();
		}
		else if (arg == "swap")
		{
			auto [v1, v2] = Pop2();
			Push({ v1, v2 });
		}
		else if (arg == "over")
		{
			auto [v1, v2] = Pop2();
			Push({ v2, v1, v2 });
		}
		else if (arg == "rot")
		{
			auto [v1, v2, v3] = Pop3();
			Push({ v2,v1,v3 });
		}
		else if (arg == "unrot")
		{
			auto [v1, v2, v3] = Pop3();
			Push({ v1,v3,v2 });
		}
		else if (arg == "depth")
		{
			Push(fmt::format("{}", nextopen));
		}
		/*
	{"roll","x1 x2.. xn n -> , like rot, but n items ",StackOp,"roll"},
	{"rolld","x1 x2 .. xn n -> , reverse of rotn",StackOp,"rolld"},
	{"pick","xn ... x1 n -> , copies item xn to top",StackOp,"pick"},
	{"unpick","X -> , NOT opposite of unpick. removes item level 1, ",StackOp,"unpick"},
		*/
		else
			throw runtime_error(fmt::format("unknown stack op {}", arg));
	}

	template<typename T>
	void Push(T item)
	{
		Resize();
		stack[nextopen++] = item;
	}
	template<typename T>
	T Pop()
	{

		auto v = Pop();
		if (!holds_alternative<T>(v))
			throw runtime_error(fmt::format("Invalid popped item type"));
		return get<T>(v);
	}
	Item Pop()
	{
		if (nextopen <= 0) throw runtime_error("empty stack");
		return stack[--nextopen];
	}

	template<typename T>
	T Peek() const
	{
		const auto & v = Peek(0);
		if (!holds_alternative<T>(v))
			throw runtime_error(fmt::format("Invalid popped item type."));
		return get<T>(stack[nextopen - 1]);
	}
	Item Peek(int n) const
	{
		auto index = nextopen - 1 - n;
		if (nextopen < 0) throw runtime_error("stack out of bounds");
		return stack[index];
	}
};