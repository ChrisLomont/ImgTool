#pragma once
#include <string>
#include "State.h"
#include "fmt/fmt/format.h"
using namespace std; // todo remove
/*------------------------ Eval and operations -------------------------*/
template<typename R> using NonaryFunc = function<R()>;
template<typename R, typename T1> using UnaryFunc = function<R(T1)>;
template<typename R, typename T1, typename T2> using BinaryFunc = function<R(T1, T2)>;
template<typename R, typename T1, typename T2, typename T3> using TrinaryFunc = function<R(T1, T2, T3)>;

/*


*/

using FuncD = NonaryFunc<double>;
using FuncDD = UnaryFunc<double, double>;
using FuncDDD = BinaryFunc<double, double, double>;
using FuncDDDD = TrinaryFunc<double, double, double, double>;
using FuncDSS = BinaryFunc<double, string, string>;
using FuncSSS = BinaryFunc<string, string, string>;

// todo - make these more table driven, Julia multi dispatch style

void Do1(State& s, const FuncDD& f)
{
	const auto v = s.Pop();
	if (holds_alternative<double>(v))
	{
		s.Push(f(get<double>(v)));
	}
	else if (holds_alternative<ImagePtr>(v))
	{
		const auto img = get<ImagePtr>(v);
		img->Apply([&](Color& c) {
			c.ApplyRGB(f);
			}
		);
		s.Push(img);
	}
	else
		throw runtime_error("Invalid types in op");
}
void Do2(State& s, const FuncDDD& f)
{
	const auto item2 = s.Pop(); // NOTE 1 and 2 swap here!
	const auto item1 = s.Pop();

	if (holds_alternative<double>(item1) && holds_alternative<double>(item2))
	{
		const auto v1 = get<double>(item1);
		const auto v2 = get<double>(item2);
		const auto ans = f(v1, v2);
		s.Push(ans);
	}
	else if (holds_alternative<double>(item1) && holds_alternative<ImagePtr>(item2))
	{
		const auto v = get<double>(item1);
		const auto src = get<ImagePtr>(item2);
		auto [w, h] = src->Size();
		const ImagePtr dst = Image::Make(w, h);
		dst->Apply([&](int i, int j) {
			auto c = src->Get(i, j);
			c.r = f(v, c.r);
			c.g = f(v, c.g);
			c.b = f(v, c.b);
			c.a = 1.0; // todo - blend?
			return c;
			}
		);
		s.Push(dst);
	}
	else if (holds_alternative<ImagePtr>(item1) && holds_alternative<double>(item2))
	{
		const auto src = get<ImagePtr>(item1);
		const auto v = get<double>(item2);
		auto [w, h] = src->Size();
		const ImagePtr dst = Image::Make(w, h);
		dst->Apply([&](int i, int j) {
			auto c = src->Get(i, j);
			c.r = f(c.r, v);
			c.g = f(c.g, v);
			c.b = f(c.b, v);
			c.a = 1.0; // todo - blend?
			return c;
			}
		);
		s.Push(dst);
	}
	else if (holds_alternative<ImagePtr>(item1) && holds_alternative<ImagePtr>(item2))
	{
		const auto img1 = get<ImagePtr>(item1);
		const auto img2 = get<ImagePtr>(item2);
		auto [w, h] = img1->Size();
		const ImagePtr img = Image::Make(w, h);
		img->Apply([&](int i, int j) {
				const auto c1 = img1->Get(i, j);
				const auto c2 = img2->Get(i, j);
			Color c;
			c.r = f(c1.r, c2.r);
			c.g = f(c1.g, c2.g);
			c.b = f(c1.b, c2.b);
			c.a = 1.0; // todo - blend?
			return c;
			}
		);
		s.Push(img);
	}
	else
		throw runtime_error("Invalid types in op");
}
void Do3(State& s, const FuncDDDD& f)
{
	const auto c = s.Pop<double>();
	const auto b = s.Pop<double>();
	const auto item = s.Pop();
	if (holds_alternative<double>(item))
	{
		const auto a = get<double>(item);
		const auto v = f(a, b, c);
		s.Push(v);
	}
	else if (holds_alternative<ImagePtr>(item))
	{
		const auto img = get<ImagePtr>(item);
		img->Apply([&](Color& co) {co.ApplyRGB([&](double v) {return f(v, b, c); }); });
		s.Push(img);
	}
	else
		throw runtime_error("Invalid types in op");
}

// types of operation functions
using funcT = variant<
	FuncD,
	FuncDD,
	FuncDDD,
	FuncDSS,
	FuncSSS,
	FuncDDDD
>;
struct OpDef
{
	string name;
	vector<funcT> funcs;
};

/*
items:
==, other compares
int int
dbl dbl
int dbl or dbl int: conv to double, ==
string string
img img

+: int int, int dbl or rev, string string

int,dbl, img math:
abs, clamp, 

dbl, img

*/

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
} // https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
const vector<OpDef> opDefs = {
	{"abs",{
		[](double v) { return abs(v); },
		//[](int v) { return abs(v); }
	}},
	{"sin",{
		[](double v) { return sin(v); },
	}},
	{"cos",{
		[](double v) { return cos(v); },
	}},
	{"==", {
		[](double a, double b) { return static_cast<double>(a == b); },
		[](string a, string b) { return static_cast<double>(a == b); },
	}},
	{"!=", {
		[](double a, double b) { return static_cast<double>(a != b); },
		[](string a, string b) { return static_cast<double>(a != b); },
	}},
	{"<=", {
		[](double a, double b) { return static_cast<double>(a <= b); },
		[](string a, string b) { return static_cast<double>(a <= b); },
	}},
	{">=", {
		[](double a, double b) { return static_cast<double>(a >= b); },
		[](string a, string b) { return static_cast<double>(a >= b); },
	}},
	{"<", {
		[](double a, double b) { return static_cast<double>(a < b); },
		[](string a, string b) { return static_cast<double>(a < b); },
	}},
	{">", {
		[](double a, double b) { return static_cast<double>(a > b); },
		[](string a, string b) { return static_cast<double>(a > b); },
	}},

	{"+", {
		[](double a, double b) { return a + b; },
		[](string a, string b) { return b + a; },
	}},
	{"-", {
		[](double a, double b) { return a - b; },
	}},
	{"*", {
		[](double a, double b) { return a * b; },
	}},
	{"/", {
		[](double a, double b) { return a / b; },
	}},
	{"mod", {
		[](double a, double b) { return fmod(a,b); }, // todo - make positive mod version?
	}},

	{"clamp",{[](double a, double b, double c) { return clamp(a,b,c); }}},

	{"round",{[](double v) { return round(v); }}},
	{"floor",{[](double v) { return floor(v); }}},
	{"ceil",{[](double v) { return ceil(v); }}},

	{"min",{[](double a, double b) { return min(a,b); }}},
	{"max",{[](double a, double b) { return max(a,b); }}},

	{"pow",{[](double a, double b) { return pow(a,b); }}},
	{"sqrt",{
		[](double v) { return sqrt(v); },
	}},

	{"exp",{[](double x) { return exp(x); }}},
	{"log",{[](double val, double base) { return log(val) / log(base); }}},

	{"neg",{[](double v) { return -v; }}},
	{"sign",{[](double v) { return sgn(v); }}},

	{"pi",{[]() { return static_cast<double>(std::numbers::pi); }}},
	{"e",{[]() { return  static_cast<double>(std::numbers::e);  }}},
};

void MathOp(State& s, const string& args)
{
	for (auto& op : opDefs)
	{
		if (op.name == args)
		{
			for (auto& vf : op.funcs)
			{
				if (holds_alternative<FuncD>(vf))
				{
					const auto f = get<FuncD>(vf);
					s.Push(f());
					return;
				}
				else if (holds_alternative<FuncDD>(vf) && s.NextTypes("X"))
				{
					Do1(s, get<FuncDD>(vf));
					return;
				}
				else if (holds_alternative<FuncDDD>(vf) && s.NextTypes("XX")) // todo- make these in op def? DD or XX
				{
					Do2(s, get<FuncDDD>(vf));
					return;
				}
				else if (holds_alternative<FuncDDDD>(vf) && s.NextTypes("XXX"))
				{
					Do3(s, get<FuncDDDD>(vf));
					return;
				}
				else if (holds_alternative<FuncSSS>(vf) && s.NextTypes("SS"))
				{
					const auto f = get<FuncSSS>(vf);
					const auto v1 = s.Pop<string>();
					const auto v2 = s.Pop<string>();
					s.Push(f(v1, v2));
					return;
				}
				else if (holds_alternative<FuncDSS>(vf) && s.NextTypes("SS"))
				{
					const auto f = get<FuncDSS>(vf);
					const auto v1 = s.Pop<string>();
					const auto v2 = s.Pop<string>();
					s.Push(f(v1, v2));
					return;
				}
#if 0
				if (holds_alternative<NonaryFunc>(v))
					Do0(s, get<NonaryFunc>(v));
				else if (holds_alternative<UnaryFunc>(v))
					Do1(s, get<UnaryFunc>(v));
				else if (holds_alternative<UnaryFunc>(v))
					Do1(s, get<UnaryFunc>(v));
				else if (holds_alternative<BinaryFunc>(v))
					Do2(s, get<BinaryFunc>(v));
				else if (holds_alternative<TrinaryFunc>(v))
					Do3(s, get<TrinaryFunc>(v));
#endif
			}
		}
	}

	throw runtime_error(fmt::format("Unknown math op {}", args));
}


void LogicOp(State& s, const string& args)
{
	//{"and", "a b -> (a and b), bitwise 'and' on integers", LogicOp},
	//{ "or","a b -> (a or b), bitwise 'or' on integers",LogicOp },
	//{ "xor","a b -> (a xor b), bitwise 'xor' on integers",LogicOp },
	//{ "not","a b -> (a not b), treating 0 as false, != 0 as true, boolean not",LogicOp },
	if (args == "and")
	{
		const auto b = s.PopInt();
		const auto a = s.PopInt();
		s.Push(a & b);
	}
	else if (args == "or")
	{
		const auto b = s.PopInt();
		const auto a = s.PopInt();
		s.Push(a | b);

	}
	else if (args == "xor")
	{
		const auto b = s.PopInt();
		const auto a = s.PopInt();
		s.Push(a ^ b);

	}
	else if (args == "not")
	{
		const auto a = s.PopInt() == 0;
		s.Push(a ? 1 : 0);
	}
	else
		throw runtime_error("Unknown op");

}