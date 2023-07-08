// simple image tool
// chris lomont 2020-2023
// resize, error metrics, gamma stuff
#include <iostream>
#include <string>
#include <cstring>
#include <filesystem>
#include <numbers>
#include <fstream>
#include <regex>

// some compilers puke, even with c++ 20, so...
//#include <format>
#include "fmt/fmt/format.h"

#include "Utils.h"
#include "Image.h"
#include "Command.h"
#include "Stack.h"
#include "State.h"
#include "Timer.h"
#include "ImageOps.h"
#include "Colorspace.h"


namespace fs = std::filesystem;
using namespace std;

// NOTES:
/*syntax:
* works on stack machine - can put images on stack, do ops
  
  TODO
    - clean command descriptions, order better
	- use double instead of string for numerical stuff on stack, less conversions
	- abstract out Do0 - Do3, abstract handlers nicer
    - string ==, != > < >= <= 
    - array of items (incl other arrays)
	- get and line args as "n getarg" to get nth arg, or maybe use rcl and special name
	- eval a string as a set of commands
	- check string construction
	- printing of endlines?
    - make CSV, would be good for running image tests
    - test script to test all, catch regressions
  X - better parser, handles comments, strings inline
  X - get files matching pattern
  X - lanczos, 
    - lanczos radial
  X - test scripts
  X - run script command
  X - run external command
    - some simple drawing - box, line, text ;)
  X - add verbosity output
    - add more filters, check
    - rotations, 
  X - crop, 
  X - expand (add border, etc), 
    - shift image around ops?
    - blit and composite using alphas?
    - basic drawing, text
  X - pixel get/set
	- pixel functions
	- more pixel ops
    - Gaussian (good fast approx, Wells, PAMI, Mar 1986), edge stuff?, list filter to apply convolution?
  X - spawn command using std::system() calls
	- bilateral filter, noise, median filter, edge stuff?
* And now we recreated Image Magic :)
*/

/*
Min viable product (x = done):
Foreach file in dir (possible recursive?)
x   open image
x  get size
x   math on size
   split into path, filename base, extension
x   crop/trim/resize whatever
x   perform image metrics (ssim, psnr, maybe all, maybe diff, abs, scale, max of cells?)
   format string with result
   append result to some internal var
   format new filename (path, name+??, extension)
x   save result(1 or more images)
save/dump total psnr, results, etc....

*/

/*-----------------------------------------------------------------------*/


void StateOp(State& s, const string& args)
{
	s.StateOp(args);
}

void StackOp(State& s, const string& args)
{
	if (s.verbosity >= 1)
		cout << fmt::format("Stack operation: {}\n", args);

	s.StackOp(args);
}

void ImageSize(State& s, const string& args)
{
	auto img = s.Peek<ImagePtr>();
	auto [w, h] = img->Size();
	s.Push(fmt::format("{}", w));
	s.Push(fmt::format("{}", h));
}

/*------------------------ Eval and operations -------------------------*/
using NonaryFunc = function<double()>;
using UnaryFunc = function<double(double)>;
using BinaryFunc = function<double(double, double)>;
using TrinaryFunc = function<double(double, double, double)>;
// todo - make these more table driven, Julia multi dispatch style
void Do0(State& s, const NonaryFunc & f)
{
	s.Push(format("{}",f()));
}
void Do1(State& s,const UnaryFunc & f)
{
	auto v = s.Pop();
	if (holds_alternative<string>(v))
	{
		s.Push(fmt::format("{}", f(ParseDouble(get<string>(v)))));
	}
	else if (holds_alternative<ImagePtr>(v))
	{
		auto img = get<ImagePtr>(v);
		img->Apply([&](Color & c) {
			c.ApplyRGB(f);
			}
		);
		s.Push(img);
	}
	else
		throw runtime_error("Invalid types in op");
}
void Do2(State& s, const BinaryFunc & f)
{
	auto item2 = s.Pop(); // NOTE 1 and 2 swap here!
	auto item1 = s.Pop();

	if (holds_alternative<string>(item1) && holds_alternative<string>(item2))
	{
		auto v1 = ParseDouble(get<string>(item1));
		auto v2 = ParseDouble(get<string>(item2));
		auto ans = f(v1, v2);
		s.Push(fmt::format("{}", ans));
	}
	else if (holds_alternative<string>(item1) && holds_alternative<ImagePtr>(item2))
	{
		auto v = ParseDouble(get<string>(item1));
		auto src = get<ImagePtr>(item2);
		auto [w, h] = src->Size();
		ImagePtr dst = make_shared<Image>(w, h);
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
	else if (holds_alternative<ImagePtr>(item1) && holds_alternative<string>(item2))
	{
		auto src = get<ImagePtr>(item1);
		auto v = ParseDouble(get<string>(item2));
		auto [w, h] = src->Size();
		ImagePtr dst = make_shared<Image>(w, h);
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
		auto img1 = get<ImagePtr>(item1);
		auto img2 = get<ImagePtr>(item2);
		auto [w, h] = img1->Size();
		ImagePtr img = make_shared<Image>(w, h);
		img->Apply([&](int i, int j) {
			auto c1 = img1->Get(i, j);
			auto c2 = img2->Get(i, j);
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
void Do3(State& s, const TrinaryFunc & f)
{
	auto c = ParseDouble(s.Pop<string>());
	auto b = ParseDouble(s.Pop<string>());
	auto item = s.Pop();
	if (holds_alternative<string>(item))
	{
		auto a = ParseDouble(get<string>(item));
		auto v = f(a, b, c);
		s.Push(fmt::format("{}", v));
	}
	else if (holds_alternative<ImagePtr>(item))
	{
		auto img = get<ImagePtr>(item);
		img->Apply([&](Color& co) {co.ApplyRGB([&](double v) {return f(v, b, c); }); });
		s.Push(img);
	}
	else
		throw runtime_error("Invalid types in op");
}


struct OpDef
{
	string name;
	variant<
		NonaryFunc,
		UnaryFunc,
		BinaryFunc,
		TrinaryFunc
	> func;
};
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
} // https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
const vector<OpDef> opDefs = {
	{"abs",[](double v) { return abs(v); }},
	{"sin",[](double v) { return sin(v); }},
	{"cos",[](double v) { return cos(v); }},
	{"+",[](double a, double b) { return a + b; }},
	{"-",[](double a, double b) { return a - b; }},
	{"*",[](double a, double b) { return a * b; }},
	{"/",[](double a, double b) { return a / b; }},
	{"mod",[](double a, double b) { return fmod(a,b); }}, // todo - make positive mod version?
	
	{"==",[](double a, double b) { return a == b; }},
	{"!=",[](double a, double b) { return a != b; }},
	{"<=",[](double a, double b) { return a <= b; }},
	{">=",[](double a, double b) { return a >= b; }},
	{"<",[](double a, double b) { return a < b; }},
	{">",[](double a, double b) { return a > b; }},
	{"clamp",[](double a, double b, double c) { return clamp(a,b,c); }},

	{"round",[](double v) { return round(v); }},
	{"floor",[](double v) { return floor(v); }},
	{"ceil",[](double v) { return ceil(v); }},
	
	{"min",[](double a, double b) { return min(a,b); }},
	{"max",[](double a, double b) { return max(a,b); }},

	{"pow",[](double a, double b) { return pow(a,b); }},
	{"exp",[](double x) { return exp(x); }},
	{"log",[](double val, double base) { return log(val)/log(base); }},

	{"neg",[](double v) { return -v; }},
	{"sign",[](double v) { return sgn(v); }},

	{"pi",[]() { return (double)std::numbers::pi; }},
	{"e",[]() { return  (double)std::numbers::e;  }},
};

// try string add, commit it if so, return true, else false
bool StringOp(State & s, const string & args)
{
	if (s.Pos() < 2) return false;
	auto i1 = s.Pop();
	auto i2 = s.Pop();
	if (holds_alternative<string>(i1) && holds_alternative<string>(i2))
	{
		auto s1 = get<string>(i1);
		auto s2 = get<string>(i2);
		if (!IsDouble(s1) && !IsDouble(s2))
		{
			if (args == "+")
			{
				s.Push(s2 + s1);
				return true;
			}
			/*
			else if (args == "==")
			{
				s.Push(s2 == s1 ? "1" : "0");
				return true;
			}
			else if (args == "!=")
			{
				s.Push(s2 != s1 ? "1" : "0");
				return true;
			}
			else if (args == "<")
			{
				s.Push(s2 < s1 ? "1" : "0");
				return true;
			}
			else if (args == ">")
			{
				s.Push(s2 > s1 ? "1" : "0");
				return true;
			}
			else if (args == "<=")
			{
				s.Push(s2 <= s1 ? "1" : "0");
				return true;
			}
			else if (args == ">=")
			{
				s.Push(s2 >= s1 ? "1" : "0");
				return true;
			}
			*/
		}
	}
	// put back
	s.Push(i2);
	s.Push(i1);
	return false;
}


void MathOp(State& s, const string& args)
{
	// special case: string ops: +, ==, !=, <, >, <=, >=
	// todo - make better way to merge these
	if (StringOp(s, args))
	{
		return;
	}

	for (auto& op : opDefs)
	{
		if (op.name == args)
		{
			auto& v = op.func;
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
			else
				throw runtime_error("Invalid func type");
			return;
		}
	}
	throw runtime_error(fmt::format("Unknown math op {}",args));
}
void Print(State& s, const string& args)
{
	int n = 1;
	if (args == "printn")
		n = ParseInt(s.Pop<string>());
	if (s.verbosity >= 2)
		cout << "Printing:\n";
	for (int i = 0; i < n; ++i)
	{
		auto v = s.Pop();
		if (s.verbosity >= 2)
			cout << fmt::format("{}: ",i);
		cout << FormatItem(v,s.verbosity>=2);
	}
	cout << endl;
}

void GetFiles(State & s, const string & args)
{
	auto regtext = s.Pop<string>();
	auto path    = s.Pop<string>();
	vector<string> files;
	regex reg(regtext);

	for (const auto& entry : fs::directory_iterator(path))
	{
		if (entry.is_regular_file())
		{
			auto path = entry.path();
			auto fn = path.filename().string();
			if (regex_match(fn, reg))
				files.push_back(path.string());
		}
	}
	sort(files.begin(), files.end());
	for (auto& f : files)
	{
		s.Push(f);
	}
	s.Push(fmt::format("{}",files.size()));
}

void PixelOp(State & s, const string & args)
{
	if (args == "getpixel")
	{
		auto j   = ParseInt(s.Pop<string>());
		auto i   = ParseInt(s.Pop<string>());
		auto img = s.Pop<ImagePtr>();
		auto c = img->Get(i,j);
		s.Push(img);
		s.Push(c.r);
		s.Push(c.g);
		s.Push(c.b);
		s.Push(c.a);
	}
	else if (args == "setpixel")
	{
		// { "setpixel", "img i j r g b a -> img, writes pixel", PixelOp },
		auto a = ParseDouble(s.Pop<string>());
		auto b = ParseDouble(s.Pop<string>());
		auto g = ParseDouble(s.Pop<string>());
		auto r = ParseDouble(s.Pop<string>());
		auto j = ParseInt(s.Pop<string>());
		auto i = ParseInt(s.Pop<string>());
		auto img = s.Pop<ImagePtr>();
		Color c(r,g,b,a);
		img->Set(i, j, c);
		s.Push(img);
	}
	else
		throw runtime_error(fmt::format("Unknown pixel op {}",args));

}


vector<Command> commands = {
	{"read","filename -> image, loads image",ReadImage},
	{"write","image filename -> ,  outputs saved image",WriteImage},
	{"colorspace","image space -> image', where space=[linear|sRGB|YCbCr|RGB], does conversion",ColorTransform},
	{"error","im1 im2 errtype -> errval, prints error type mse, psnr, ssim",ImageError},
	{"maxc","img -> max, max value of all r,g,b values in image",ImageError},
	{"size","img -> w h, where w,h is size in pixels",ImageSize},
	{"files","path regex -> f1 f2 ... fn n, reads files matching regex, pushes on stack with count",GetFiles},

	{"getpixel","img i j -> img r g b a, reads pixel 0-1",PixelOp},
	{"setpixel","img i j r g b a -> img, writes pixel 0-1",PixelOp},

	// and,or,not,xor,rand, rdz(randseed), >>,<<	
	// ticks = time, 
	// type - object type


	{"abs","item -> abs(img)",MathOp},
	{"ceil","item -> ceil(item)",MathOp},
	{"floor","item -> floor(img)",MathOp},
	{"round","item -> round(item)",MathOp},
	{"min","a b -> min(a,b)",MathOp},
	{"max","a b -> max(a,b)",MathOp},
	{"clamp","item1 a b -> clamp(item1,a,b)",MathOp},

	{"sin","item -> abs(img)",MathOp},
	{"cos","item -> abs(img)",MathOp},

	{"pi"," -> pi",MathOp},
	{"e","  -> e",MathOp},

	{"pow","item1 item2 -> pow(item1,item2)",MathOp},
	{"exp","a -> e^a ",MathOp},
	{"log","val base -> log_base(val)",MathOp},

	{"neg","a -> -a",MathOp},
	{"sign","a -> sign(a), is -1,0,1",MathOp},

	{"+","item1 item2 -> item1+item2",MathOp},
	{"-","item1 item2 -> item1-item2",MathOp},
	{"*","item1 item2 -> item1*item2",MathOp},
	{"/","item1 item2 -> item1/item2",MathOp},
	{"mod","a b -> a mod b",MathOp},
	{"==","item1 item2 -> item1==item2, 0 if false, else 1",MathOp},
	{"!=","item1 item2 -> item1!=item2, 0 if false, else 1",MathOp},
	{">=","item1 item2 -> item1>=item2, 0 if false, else 1",MathOp},
	{"<=","item1 item2 -> item1<=item2, 0 if false, else 1",MathOp},
	{">","item1 item2 -> item1>item2, 0 if false, else 1",MathOp},
	{"<","item1 item2 -> item1<item2, 0 if false, else 1",MathOp},

	// image ops
	{"resize","img w h style -> img', resize to w h by style nn,bilinear,bicubic,lanczos2,lanczos3,lanczos4",ResizeImage},
	{"resize%","img v style -> img', resize by v%, style as above",ResizeImage},
	{"resize*","img m style -> img', resize by multiplier m, style as above",ResizeImage},

	{"gaussian","img s -> img' , gaussian blur, std dev s",GaussianBlur},

	{"rotate","TODO: img angle expand -> img', rotate image by angle degrees, expand true makes bigger to center, false keeps size",RotateImage},

	{"crop","img x1 y1 x2 y2 -> img', crop image to rectangle (x1,y1)-(x2,y2) inclusive", CropImage},
	{"pad", "img top bottom left right r g b a -> img2, pad image with given color, given pixel margins", PadImage},
	{"flipx", "img -> img2, flip image", FlipImage},
	{"flipy", "img -> img2, flip image", FlipImage},
	// todo - draw, text, trim


	// stack commands
	{"dup","a -> a a, duplicates top item",StackOp},
	{"dup2","a b -> a b a b, duplicates top item",StackOp},
	{"dupn","x1 x2 .. xn n -> x1 x2 .. xn x1 x2 .. xn, duplicate top n",StackOp},
	{"drop","a -> , drops top item",StackOp},
	{"drop2","a b -> , drops top item",StackOp},
	{"dropn","x1 x2 .. xn n -> , drops top n",StackOp},
	{"swap","a b -> b a , swaps top 2 items",StackOp},
	{"over","a b -> a b a , copies item at level 2 to top",StackOp},
	{"rot","3 2 1 -> 2 1 3 , rotates item in level 3 to level 1, 1 to 2, 2 to 3",StackOp},
	{"unrot","3 2 1 -> 1 3 2 , rotates opposite of rot",StackOp},
	{"roll","x1 x2.. xn n -> x2 x3 ... xn x1  , like rot, but n items ",StackOp},
	{"rolld","x1 x2 .. xn n -> xn x1 x2 x3 ... xn-1 , reverse of roll",StackOp},
	{"pick","xn ... x1 n -> xn ... x1 xn , copies item xn to top",StackOp},
// todo - ??	{"unpick","X -> , NOT opposite of unpick. removes item level 1",StackOp},
	{"depth","... -> n , pushes depth to top of stack",StackOp},

	{"print","item -> , prints top item",Print},
	{"printn","x1 x2 ... xn  n -> , prints top N items",Print},

	// flow & state
	{"label", "name -> , create named label for next item index",StateOp},
	{"ja", "label -> , JumpAlways: goto label",StateOp},
	{"je", "val label -> , JumpEqual: if val != 0, goto label",StateOp},
	{"halt" ,"-> , stops program", StateOp},
	{"sto" ,"item name -> , store item in name", StateOp},
	{"rcl" ,"name -> item, look up item", StateOp},
	{"dumpstate" ," -> , print out state items", StateOp},
	{"system" ,"cmd -> return_value, execute cmd on system call - WARNING - be careful!", StateOp},
	{"verbosity","v -> , set verbosity 0=none, 1=info, 2= all", StateOp},
	{"if","t1 t2 .. tn f1 f2 .. fm n m b -> ti or fj, if b != 0, keep t1..tn, else keep f1..fm", StateOp},

	// loop
	{"rangeloop","min max -> , loops over index in [min,max], each iter puts index on stack, use endloop",StateOp},
	{"itemloop","i1 i2 .. in n -> , loops over items in {i1,i2,..,in}, each iter puts item then index i=0+ on stack, use endloop", StateOp},
	{"endloop"," -> , ends loop, jumps to top", StateOp},

	// store items in array
	

	// str-> (execute string), 
	// ->str (object to string?), 
	// vars - dump stored items (vars, labels)
    //	if/then/else, switch?
};


void ShowUsage()
{
	cout << "Usage: This is an RPN based image tool. Command args are RPN commands.\n";
	cout << "       Commands either on command line or run as --script filename\n";
	cout << "       --verbose to print more\n";
	cout << "       Each command shows what it does to the stack.\n";
	for (auto& c : commands)
	{
		cout << fmt::format("{:12}: {}\n", c.name, c.description);
	}
}

bool Process(const vector<string> & tokens, bool verbose)
{
	State state;
	state.verbosity = verbose?2:1;
	
	state.programPosition = 0; 
	try {
		while (state.programPosition < tokens.size())
		{
			const string & token = tokens[state.programPosition];
			state.programPosition++; // next position

			int cIndex = 0;
			while (cIndex < commands.size())
			{
				auto& c = commands[cIndex];
				if (c.name == token)
				{
					if (state.verbosity >= 2)
						cout << fmt::format("Executing {}\n",c.name);
					c.op(state, token);
					break;
				}
				++cIndex;
			}
			if (cIndex >= commands.size())
			{ // was no item, push it
				state.Push(token);
			}

		}
	}
	catch (runtime_error & e)
	{
		cerr << fmt::format("Exception {}\n", e.what());
		return false;
	}
	catch (...)
	{
		cerr << fmt::format("Exception: should not reach this handler!");
		return false;
	}
	cout << fmt::format("Stack depth on exit {}:",state.Pos());
	return true;
}


void GetScriptTokens(vector<string> & tokens, const string & filename)
{
	std::ifstream file(filename);
	std::string line;
	while (std::getline(file, line))
	{
		line = trim(line);

		// string split by regex
		// grab strings "..." whole
		// take comment #... to end of line
		// split on spaces
		// ignore multiple spaces

		regex rgx("(#.+$|\"[^\"]+\"|[^\\s]+)");
		sregex_token_iterator iter(
			line.begin(),
			line.end(),
			rgx);
		sregex_token_iterator end;
		for (; iter != end; ++iter)
		{
			string token = *iter;
			//cout << format("tok: <{}>",token);
			if (token[0] == '#') continue; // comment
			if (token[0] == '"') token = token.substr(1, token.size() - 2); // string
			//cout << format(" => <{}>\n", token);
			tokens.push_back(token);
		}
	}
}

int main(int argc, char ** argv)
{
	cout << "Chris Lomont's RPN image tool v0.1\n";
	if (argc <= 1)
	{
		ShowUsage();
		return -1;
	}
	vector<string> tokens;
	bool verbose = false;
	
	int argpos = 1; // skip initial exe
	while (argpos < argc&& argv[argpos][0]=='-')
	{ // parse options
		string opt(argv[argpos]);
		argpos++;
		if (opt == "--script")
		{
			cout << "Executing script " << argv[argpos] << endl;
			GetScriptTokens(tokens, argv[argpos++]);
		}
		else if (opt == "--verbose")
		{
			cout << "Verbose = true\n";
			verbose = true;
		}
		else {
			cerr << fmt::format("Unknown option {}\n",opt);
			return -1;
		}
	}
	// tokenize any other command line options
	for (int i = argpos; i < argc; ++i)
		tokens.push_back(argv[i]);

	cout << "Current path is " << fs::current_path() << '\n'; 
	auto retval = Process(tokens, verbose) ? 1 : 0;
	cout << "Done\n";
	return retval;
}

