
// simple image tool
// chris lomont 2020-2023
// resize, error metrics, gamma stuff

// https://github.com/ChrisLomont/ImgTool
// to compile on *nix, g++ main.cpp fmt/format.cc -std=c++20 -O2 -o imgtool
// then copy imgtool to /usr/local/bin to make available on terminals

#include <iostream>
#include <string>
#include <cstring>
#include <filesystem>
#include <numbers>
#include <fstream>
#include <regex>
#include <typeinfo>
#include <set>

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
#include "MathOps.h"
#include "csv.h"

namespace fs = std::filesystem;
using namespace std;

const int VERSION_MAJOR = 0;
const int VERSION_MINOR = 2;

// NOTES:
/*syntax:
* works on stack machine - can put images on stack, do ops

    - doc assumptions: read file does no conversion except byte values -> double rgba
 
  
  TODO
    - 
    - string with func name runs, e.g., "error" triggers error. Make quoted strings allowed as strings and labels and var names...
    - better help - explin images loaded as is, more
    - point to github
      - add complete color space support - for experiments
          - tag images with color space details to allow proper "to" conversions
          - how to read/save such images? stb likely does not handle it....
  X - set edge boundary mode for image stuff: reflection, clamp, etc...
    - filename processing, or more general regex stuff?
	- abstract out rpn and stack engine, it's decent
	- script processing - look direct name, if not, and name has no /, then look in scripts subdir
	- default scripts, ex: diffimg, errimg, histogram, etc.
  X	- script can see # of args passed
    - better io control to make silent via verbosity
	- console tools run offline can hide output
  X - function support: <name> subroutine/endsub and <name> gosub/return
  X	- printing of endlines
	- format prints
    - assert
	- type on stack command
	- eval a string as a set of commands
	- type converters ,
  X    ->str,
	   ->float, 
	- replace Hoppe code, license not clear
	- clean command descriptions, order better
  X - use double instead of string for numerical stuff on stack, less conversions
	- abstract out Do0 - Do3, abstract handlers nicer
  X - string ==, != > < >= <= 
    - array of items (incl other arrays)
  X	- get and line args as "n getarg" to get nth arg, or maybe use rcl and special name
	- check string construction
  X - make CSV, would be good for running image tests
    - test script to test all, catch regressions
  X - better parser, handles comments, strings inline
  X - get files matching pattern
  X - lanczos, 
  X - lanczos radial
  X - test scripts
  X - run script command
  X - run external command
  X - add verbosity output
    - add more filters, check
    - rotations, 
  X - crop, 
  X - expand (add border, etc), 
    - shift image around ops?
  X - blit and composite using alphas?

  X	- version
  X - int rgba <-> double rgba

  X	- make new image, blank,
  X - pixel get/set
	- pixel functions
	- more pixel ops
	- draw line
	- draw rect
	- draw circle
	- text
	- fill, ops

  X - Gaussian 
	   - (good fast approx, Wells, PAMI, Mar 1986), edge stuff?, list filter to apply convolution?
  X - spawn command using std::system() calls
	- bilateral filter, noise, median filter, edge stuff?

* And now we recreated Image Magic :)
*/

/*
Min viable product (x = done):
x Foreach file in dir (possible recursive?)
x   open image
x  get size
x   math on size
   split into path, filename base, extension
x   crop/trim/resize whatever
x   perform image metrics (ssim, psnr, maybe all, maybe diff, abs, scale, max of cells?)
  x format string with result
  x append result to some internal var
  x format new filename (path, name+??, extension)
x   save result(1 or more images)
x save/dump total psnr, results, etc....

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

void RandOp(State& s, const string& args)
{
	if (args == "rand")
	{
		auto b = s.PopInt();
		auto a = s.PopInt();
		auto val = randUniform(s.randState, a, b);
		s.Push(val);
	}
	else if (args == "srand")
	{
		auto seed = s.PopInt();
		s.randState = randSeed(seed);
	}
	else throw runtime_error(fmt::format("Unknown command {}",args));
}

void Print(State& s, const string& args)
{
	if (args == "endl")
	{
		s.Push("\n");
		return;
	}
	int n = 1;
	if (args == "printn")
		n = s.PopInt();
	if (s.verbosity >= 2)
		cout << "Printing:\n";
	for (int i = 0; i < n; ++i)
	{
		auto v = s.Pop();
		if (s.verbosity >= 2)
			cout << fmt::format("{}: ",i);
		cout << FormatItem(v,s.verbosity>=2) << " ";
	}
	cout << endl;
}
void SystemOp(State& s, const string& args)
{
		if (args == "version") 
		{
			// version
			s.Push(VERSION_MAJOR);
			s.Push(VERSION_MINOR);
		}
		else if (args == "timeus")
		{
			s.Push(s.elapsedUs());
		}
		else if (args == "arg")
		{
			auto n = s.PopInt();
			s.Push(s.args[n]);
		}
		else if (args == "argcount")
		{
			s.Push((int)s.args.size());
		}
		else throw runtime_error("Unknown system op");
}

void LogicOp(State& s, const string& args)
{
	//{"and", "a b -> (a and b), bitwise 'and' on integers", LogicOp},
	//{ "or","a b -> (a or b), bitwise 'or' on integers",LogicOp },
	//{ "xor","a b -> (a xor b), bitwise 'xor' on integers",LogicOp },
	//{ "not","a b -> (a not b), treating 0 as false, != 0 as true, boolean not",LogicOp },
	if (args == "and")
	{

	}
	else if (args == "and")
	{
		auto b = s.PopInt();
		auto a = s.PopInt();
		s.Push(a & b);
	}
	else if (args == "or")
	{
		auto b = s.PopInt();
		auto a = s.PopInt();
		s.Push(a | b);

	}
	else if (args == "xor")
	{
		auto b = s.PopInt();
		auto a = s.PopInt();
		s.Push(a ^ b);

	}
	else if (args == "not")
	{
		auto a = s.PopInt() == 0;
		s.Push(a ? 1 : 0);
	}
	else 
		throw runtime_error("Unknown op");

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
	s.Push((int)files.size());
}

void PixelOp(State & s, const string & args)
{
	if (args == "getpixel")
	{
		auto j   = s.PopInt();
		auto i   = s.PopInt();
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
		auto a = s.Pop<double>();
		auto b = s.Pop<double>();
		auto g = s.Pop<double>();
		auto r = s.Pop<double>();
		auto j = s.PopInt();
		auto i = s.PopInt();
		auto img = s.Pop<ImagePtr>();
		Color c(r,g,b,a);
		img->Set(i, j, c);
		s.Push(img);
	}
	else
		throw runtime_error(fmt::format("Unknown pixel op {}",args));

}

set<string> filesRead;

void GetScriptTokens(vector<string>& tokens, const string& filename)
{
	// insert each line
	if (!fs::exists(filename))
	{
		throw runtime_error(fmt::format("File {} does not exist", filename));
	}
	if (filesRead.find(filename) != filesRead.end())
	{
		return; // do not load again
	}
	cout << "Reading file " << filename << endl;

	filesRead.insert(filename); // mark first to prevent recursive includes


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

void IncludeFile(State & s, const string & args)
{
//	{"include", " filename -> , include file as text, each file included at most once", IncludeFile},
	if (args == "include")
	{
		auto filename = s.Pop<string>();
		
		vector<string> tokens;
		GetScriptTokens(tokens, filename);
		s.tokens.insert(next(s.tokens.begin(),s.programPosition), tokens.begin(), tokens.end());
	}
	else
		throw runtime_error("unknown Include option");
}

vector<Command> commands = {
	// image stuff
	{"read","filename -> image, loads image",ImageOp},
	{"write","image filename -> ,  outputs saved image",ImageOp},
	{"image","w h r g b a -> image, makes image size w x h, color rgba in 0-1",ImageOp},
	{"getpixel","img i j -> img r g b a, reads pixel 0-1",PixelOp},
	{"setpixel","img i j r g b a -> img, writes pixel 0-1",PixelOp},

	{"colorspace","image space -> image', where space=[linear|sRGB|YCbCr|RGB], does conversion",ColorTransform},

	{"error","im1 im2 errtype -> im1 im2 errval, prints error, errtype mse, psnr, ssim",ImageError},
	{"maxc","img -> max, max value of all r,g,b values in image",ImageError},
	{"size","img -> w h, where w,h is size in pixels",ImageOp},

	// image ops
	{"resize","img w h style -> img', resize to w h by style nn,bilinear,bicubic,lanczos2,lanczos3,lanczos4,lanczos2r,lanczos3r,lanczos4r",ResizeImage},
	{"resize%","img v style -> img', resize by v%, style as above",ResizeImage},
	{"resize*","img m style -> img', resize by multiplier m, style as above",ResizeImage},

	{"gaussian","img radius -> img' , gaussian blur with given radius",GaussianBlur},

	{"rotate","TODO: img angle expand -> img', rotate image by angle degrees, expand true makes bigger to center, false keeps size",RotateImage},

	{"crop","img x1 y1 x2 y2 -> img', crop image to rectangle (x1,y1)-(x2,y2) inclusive", CropImage},
	{"pad", "img top bottom left right r g b a -> img2, pad image with given color, given pixel margins", PadImage},
	{"flipx", "img -> img2, flip image", FlipImage},
	{"flipy", "img -> img2, flip image", FlipImage},

	{"blit", "dst src -> dst', copy pixels from src to dst", ImageOp},
	{"blitc", "dst dx dy src -> dst' copy src pixels to dst, placing dest corner at dx dy", ImageOp },
	{"blitr", "dst dx dy src x1 y1 w h  -> dst', copy rect from src x1 y1 w h to dst at dx dy", ImageOp },

	{"boundary", "img mode -> img', set sample boundary mode to clamp, reflect, reverse, tile", ImageOp },


	// todo - draw, text, trim

	{"f->i","f1 f2 .. fn n -> i1 i2 .. in, converts n values in 0-1 to n values in 0-255, useful for colors",ImageOp},
	{"i->f","i1 i2 .. in n -> f1 f2 .. fn, converts n values in 0-255 to n values in 0-1, useful for colors",ImageOp},
	//{"apply","img funcname -> img', applies function funcname(i,j,r,g,b,a)->(r,g,b,a) to image pixels.",ImageOp},

	// system
	{"files","path regex -> f1 f2 ... fn n, reads files matching regex, pushes on stack with count",GetFiles},
	{"version"," -> major minor, get version", SystemOp},
	{"timeus"," -> time_us, get elapsed time in us", SystemOp},
	{"rand","a b -> rand32(a,b), uniform random integer in [a,b)",RandOp},
	{"srand","seed -> , set random seed to integer seed",RandOp},
	{"arg", " n -> arg, get command line arg n, passed via -a item, n = 1,2,...",SystemOp},
	{"argcount", "  -> argcount, count of command line args passed via via -a",SystemOp},
	
	{"include", " filename -> , include file as text, each file included at most once",IncludeFile},

	// >>,<<	
	// type - object type

	// CSV 
	{"csvstart", " csvname header1 header2 ... n -> , start a CSV file with given headers",CsvOp},
	{"csvput"  , " val header csvname -> , stores val under header name in named csv",CsvOp},
	{"csvwrite", " csvname filename -> , ",CsvOp},

	// math
	{"abs","item -> abs(img)",MathOp},
	{"ceil","item -> ceil(item)",MathOp},
	{"floor","item -> floor(img)",MathOp},
	{"round","item -> round(item)",MathOp},
	{"min","a b -> min(a,b)",MathOp},
	{"max","a b -> max(a,b)",MathOp},
	{"clamp","item a b -> clamp(item,a,b)",MathOp},
	{"sin","item -> sin(item), values in radians",MathOp},
	{"cos","item -> cos(item), values in radians",MathOp},
	{"pi"," -> pi",MathOp},
	{"e","  -> e",MathOp},
	{"pow","item1 item2 -> pow(item1,item2)",MathOp},
	{"exp","item -> e^item ",MathOp},
	{"log","val base -> log_base(val)",MathOp},
	{"neg","item -> -item",MathOp},
	{"sign","item -> sign(item), is -1,0,1",MathOp},
	{"+","item1 item2 -> item1 + item2",MathOp},
	{"-","item1 item2 -> item1 - item2",MathOp},
	{"*","item1 item2 -> item1 * item2",MathOp},
	{"/","item1 item2 -> item1 / item2",MathOp},
	{"mod","item1 item2 -> item1 mod item2",MathOp},
	{"==","item1 item2 -> item1 == item2, 0 if false, else 1",MathOp},
	{"!=","item1 item2 -> item1 != item2, 0 if false, else 1",MathOp},
	{">=","item1 item2 -> item1 >= item2, 0 if false, else 1",MathOp},
	{"<=","item1 item2 -> item1 <= item2, 0 if false, else 1",MathOp},
	{">","item1 item2 -> item1 > item2, 0 if false, else 1",MathOp},
	{"<","item1 item2 -> item1 < item2, 0 if false, else 1",MathOp},
	{"and","a b -> (a and b), bitwise 'and' on integers",LogicOp},
	{"or","a b -> (a or b), bitwise 'or' on integers",LogicOp},
	{"xor","a b -> (a xor b), bitwise 'xor' on integers",LogicOp},
	{"not","a -> (not a), treating 0 as false, != 0 as true, boolean not",LogicOp},

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

	// printing
	{"print","item -> , prints top item",Print},
	{"printn","x1 x2 ... xn  n -> , prints top N items",Print},
	{"endl", "-> endline, pushes an endline string", Print},

	// flow & state
	{"label", "name -> , create named label for next item index",StateOp},
	{"ja", "label -> , JumpAlways: goto label",StateOp},
	{"je", "item1 item2 label -> , jump to label if item1 == item2",StateOp},

	{ "jne", "item1 item2 label -> , jump to label if item1 != item2",StateOp },
	{ "jlt", "item1 item2 label -> , jump to label if item1 < item2",StateOp },
	{ "jle", "item1 item2 label -> , jump to label if item1 <= item2",StateOp },
	{ "jgt", "item1 item2 label -> , jump to label if item1 > item2",StateOp },
	{ "jge", "item1 item2 label -> , jump to label if item1 >= item2",StateOp },

	{"halt" ,"exitcode -> , stops program, returns code", StateOp},
	{"sto" ,"item name -> , store item in name", StateOp},
	{"rcl" ,"name -> item, look up item", StateOp},
	{"dumpstate" ," -> , print out state items", StateOp},
	{"system" ,"cmd -> return_value, execute cmd on system call - WARNING - be careful!", StateOp},
	{"verbosity","v -> , set verbosity 0=none, 1=info, 2=all", StateOp},
	{"if","t1 t2 .. tn f1 f2 .. fm n m b -> ti or fj, if b != 0, keep t1..tn, else keep f1..fm", StateOp},
	{"->str","item -> 'item', formats item as string", StateOp},


	// loop
	{"rangeloop","min max -> , loops over index in [min,max], each iter puts index on stack, use endloop",StateOp},
	{"itemloop","i1 i2 .. in n -> , loops over items in {i1,i2,..,in}, each iter puts item then index i=0+ on stack, use endloop", StateOp},
	{"endloop"," -> , ends loop, jumps to top", StateOp},

	// subroutines
	{ "subroutine"," name -> , starts subroutine, ends with endsub", StateOp },
	{ "endsub"," name -> , ends subroutine", StateOp },
	{ "gosub"," name -> , jumps to subroutine ", StateOp },
	{ "return"," -> , returns from subroutine", StateOp },



	// store items in array
	

	// str-> (execute string), 
	// ->str (object to string?), 
	// vars - dump stored items (vars, labels)
    //	if/then/else, switch?
};


void ShowUsage()
{
	cout << "Usage: This is an RPN based image tool. Command args are RPN commands.\n";
	cout << "       Commands either on command line or run as -s filename\n";
	cout << "       --verbose to print more, 0=none, 1=info, 2=all\n";
	cout << "       Each command shows what it does to the stack.\n";
	for (auto& c : commands)
	{
		cout << fmt::format("{:12}: {}\n", c.name, c.description);
	}
}

Item ToItem(const string& text)
{
	if (IsDouble(text))
		return Item(ParseDouble(text));
	return Item(text);
}

bool Process(State & state, bool verbose)
{
	state.verbosity = verbose?2:1;
	
	state.programPosition = 0; 
	try {
		while (state.programPosition < state.tokens.size())
		{
			const string & token = state.tokens[state.programPosition];
			state.programPosition++; // next position

			if (state.inSubroutineDefinition && token != "endsub")
				continue; // do nothing

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
				state.Push(ToItem(token));
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

int main(int argc, char ** argv)
{
	cout << fmt::format("Chris Lomont's RPN image tool v{}.{}\n", VERSION_MAJOR, VERSION_MINOR);
	if (argc <= 1)
	{
		ShowUsage();
		return -1;
	}
	bool verbose = false;
	
	State s;

	int argpos = 1; // skip initial exe
	while (argpos < argc&& argv[argpos][0]=='-')
	{ // parse options
		string opt(argv[argpos]);
		argpos++;
		if (opt == "-s")
		{
			cout << "Executing script " << argv[argpos] << endl;
			GetScriptTokens(s.tokens, argv[argpos++]);
		}
		else if (opt == "--verbose")
		{
			cout << "Verbose = true\n";
			verbose = true;
		}
		else if (opt == "-a")
		{
			string t(argv[argpos++]);
			s.args.push_back(ToItem(t));
		}
		else {
			cerr << fmt::format("Unknown option {}\n",opt);
			return -1;
		}
	}
	// tokenize any other command line options
	for (int i = argpos; i < argc; ++i)
		s.tokens.push_back(argv[i]);

	cout << "Current path is " << fs::current_path() << '\n'; 
	auto retval = Process(s, verbose) ? 1 : 0;
	cout << "Done\n";
	return retval;
}

