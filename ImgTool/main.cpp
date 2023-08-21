
// simple image tool
// chris lomont 2020-2023
// resize, error metrics, gamma stuff

// https://github.com/ChrisLomont/ImgTool
// to compile on *nix, g++ main.cpp fmt/format.cc -std=c++20 -O2 -o imgtool
// then copy imgtool to /usr/local/bin to make available on terminals

#include <iostream>
#include <string>

// some compilers puke (clang!), even with c++ 20, so replace <format> with older one
//#include <format>
#include "fmt/fmt/format.h"

//#include <typeinfo>
//#include "Utils.h"
#include "Stack.h"
//#include "Command.h"
//#include "State.h"

#include "RPNLang.h"

#include "ImgLibInterface.h"
#include "MathOps.h"
#include "csv.h"

using namespace std;

const int VERSION_MAJOR = 0;
const int VERSION_MINOR = 3;



// NOTES:
/*syntax:
* works on stack machine - can put images on stack, do ops

    - document assumptions:
      - read file does no conversion except byte values -> double rgba
      - user must do gamma<->linear, etc.
    
 
  
  TODO
    - add features to set image format stuff: color space, alpha, gamma, etc. RGB, YUV? YCbCr?    
    - pre-multiplied alpha and correct compositing and tests
    - on exception, publish stack, token neighborhood in program, prog line/file?
    - add more ifs? if_gosub, ife_gosub (if = if {}, ife = if/else), maybe ifn, switch, etc?
  X - replace all std::Image::Make with Image::Make items
    - all ImgLib to namespace Lomont::ImgLib
    - all ImgLib to image_exception
    - pixel zoom: each pixel becomes NxN for showing images
    - gif or other animated form would be nice
    - thick lines on draw?
    - crop to bounds on image
    - ensure shift works on fractional pixels, do test against results in papers
	- do major refactor:
	 X * split out RPN engine
	     * make nicer interface for +,-, etc between RPN lang and Image interop
	     * make item type extensible between RPN and Image addition
	   * make commands infra cleaner
	 X * split out all graphics into standalone lib
	   * move many files to cpp, integrate lots of my old code
    - add virtual slice, make filters use them, unify filter stuff, template them all? :)
    - add pixel tests: idea is to make some test values, roundtrip all the things (alpha, belnding, ??)
    -  check blending works for all drawing functions
	   - do gamma tests, show errors under different things
	- {"apply","img funcname -> img', applies function funcname(i,j,r,g,b,a)->(r,g,b,a) to image pixels.",ImageOp},
	- todo - trim to pixels not matching some value
    - string with func name runs, e.g., "error" triggers error. Make quoted strings allowed as strings and labels and var names...
    - better help - explain images loaded as is, more
  X - point to github
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
	- type converters 
  X    ->str,
	   ->float, 
	- replace Hoppe code, license not clear
	- clean command descriptions, order better
  X - use double instead of string for numerical stuff on stack, less conversions
	- abstract out Do0 - Do3, abstract handlers nicer
  X - string ==, != > < >= <= 
    - array of items (incl other arrays)
       - call a list, like in HP RPN, { starts array, } closes,
       - add sort, length, append, insert, read nth
       - add items, n ->list makes a list
       - add list-> pops items out, n on top
	   - make loop work over list
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
  X - crop, 
  X - expand (add border, etc), 
  X - shift image around ops?
  X	- rotations,
  X - blit and composite using alphas?

  X	- version
  X - int rgba <-> double rgba

  X	- make new image, blank,
  X - pixel get/set
	- pixel functions
	- more pixel ops
  X - draw line
  X - draw rect
  X - draw circle
  X - text
  X - fill, 
	- ops

  X - Gaussian 
	   - (good fast approx, Wells, PAMI, Mar 1986), edge stuff?, list filter to apply convolution?
  X - spawn command using std::system() calls
	- bilateral filter, noise, median filter, edge stuff?


// some more commands?
	// >>,<<
	// type - object type

	// store items in array
	// str-> (execute string),
	// ->str (object to string?),
	// vars - dump stored items (vars, labels)
	// switch?


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




// todo - move command defs elsewhere, into respective locations
vector<Command> imageCommands = {
	// image stuff
	{"read","filename -> image, loads image",ImageOp},
	{"write","image filename -> ,  outputs saved image",ImageOp},
	{"image","w h r g b a -> image, makes image size w x h, color rgba in 0-1",ImageOp},
	{"getpixel","img i j -> img r g b a, reads pixel 0-1",ImageOp},
	{"setpixel","img i j r g b a -> img, writes pixel 0-1",ImageOp},

	{"colorspace","image space -> image', where space=[linear|sRGB|YCbCr|RGB], does conversion",ImageOp},
	{"alpha*","image -> image', applies alpha pre-multiplication",ImageOp},
	{"alpha/","image -> image', reverses alpha pre-multiplication (0 alpha -> 0,0,0,0 color)",ImageOp},

	{"error","im1 im2 errtype -> im1 im2 errval, prints error, errtype mse, psnr, ssim",ImageOp},
	{"maxc","img -> max, max value of all r,g,b values in image",ImageOp},
	{"size","img -> w h, where w,h is size in pixels",ImageOp},

	// image ops
	{"resize","img w h style -> img', resize to w h by style nn,bilinear,bicubic,lanczos2,lanczos3,lanczos4,lanczos2r,lanczos3r,lanczos4r",ImageOp},
	{"resize%","img v style -> img', resize by v%, style as above",ImageOp},
	{"resize*","img m style -> img', resize by multiplier m, style as above",ImageOp},

	{"gaussian","img radius -> img' , gaussian blur with given radius",ImageOp},

	{"rotate","img angle filter -> img', rotate image by angle degrees using filter nn,bilinear,bicubic",ImageOp},
	{"shift","img dx dy filter -> img', shift image by dx dy using filter (todo all nn for now)",ImageOp},

	{"crop","img x1 y1 x2 y2 -> img', crop image to rectangle (x1,y1)-(x2,y2) inclusive", ImageOp},
	{"pad", "img top bottom left right r g b a -> img2, pad image with given color, given pixel margins", ImageOp},
	{"flipx", "img -> img2, flip image", ImageOp},
	{"flipy", "img -> img2, flip image", ImageOp},

	{"blit", "src dst -> dst', copy pixels from src to dst", ImageOp},
	{"blitc", "src dst dx dy -> dst' copy src pixels to dst, placing dest corner at dx dy", ImageOp },
	{"blitr", "src x1 y1 w h dst dx dy -> dst', copy rect from src x1 y1 w h to dst at dx dy", ImageOp },
	{"blitover", "src dst -> dst dx dy', alpha blend src OVER dst, at dx dy", ImageOp },

	{"boundary", "img [r g b a] mode -> img', set sample boundary mode to color (with rgba), clamp, reflect, reverse, tile", ImageOp },

	{"f->i","f1 f2 .. fn n -> i1 i2 .. in, converts n values in 0-1 to n values in 0-255, useful for colors",ImageOp},
	{"i->f","i1 i2 .. in n -> f1 f2 .. fn, converts n values in 0-255 to n values in 0-1, useful for colors",ImageOp},
};

vector<Command> drawCommands = {
	{"line",    "img x1 y1 x2 y2 r g b a -> img with line",DrawOp},
	{"circle",  "img x1 y1 radius r g b a -> img with circle",DrawOp},
	{"circlef", "img x1 y1 radius r g b a -> img with filled circle",DrawOp},
	{"rect",    "img x1 y1 x2 y2 r g b a -> img with rectangle",DrawOp},
	{"rectf",   "img x1 y1 x2 y2 r g b a -> img with filled rectangle",DrawOp},
	{"text",    "img x1 y1 r g b a text 0 m -> img x2 y2, draws text in font (always 0), pixel size m, returns img and final position",DrawOp },
};



vector<Command> csvCommands = {
	// CSV 
	{"csvstart", " csvname header1 header2 ... n -> , start a CSV file with given headers",CsvOp},
	{"csvput"  , " val header csvname -> , stores val under header name in named csv",CsvOp},
	{"csvwrite", " csvname filename -> , ",CsvOp},
};


vector<Command> mathCommands = {
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
	{"sqrt","item1 -> sqrt(item)",MathOp},
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
};

vector<Command> logicCommands = {
	{"and","a b -> (a and b), bitwise 'and' on integers",LogicOp},
	{"or","a b -> (a or b), bitwise 'or' on integers",LogicOp},
	{"xor","a b -> (a xor b), bitwise 'xor' on integers",LogicOp},
	{"not","a -> (not a), treating 0 as false, != 0 as true, boolean not",LogicOp},
};



void ShowUsage(const RPNLanguage & rpnProcessor)
{
	cout << "Usage: This is an RPN based image tool. Command args are RPN commands.\n";
	cout << "       Commands either on command line or run as -s filename\n";
	cout << "       --verbose to print more, 0=none, 1=info, 2=all\n";
	cout << "       Each command shows what it does to the stack.\n";
	for (auto& c : rpnProcessor.commandList)
	{
		cout << fmt::format("{:12}: {}\n", c.name, c.description);
	}
}


int main(int argc, char** argv)
{
	cout << fmt::format("Chris Lomont's RPN image tool v{}.{}, https://github.com/ChrisLomont/ImgTool\n", VERSION_MAJOR, VERSION_MINOR);

	RPNLanguage rpnProcessor(VERSION_MAJOR, VERSION_MINOR);

	// prepare commands	
	rpnProcessor.commandList.insert(rpnProcessor.commandList.end(), imageCommands.begin(), imageCommands.end());
	rpnProcessor.commandList.insert(rpnProcessor.commandList.end(), drawCommands.begin(), drawCommands.end());
	rpnProcessor.commandList.insert(rpnProcessor.commandList.end(), csvCommands.begin(), csvCommands.end());
	rpnProcessor.commandList.insert(rpnProcessor.commandList.end(), mathCommands.begin(), mathCommands.end());
	rpnProcessor.commandList.insert(rpnProcessor.commandList.end(), logicCommands.begin(), logicCommands.end());
	rpnProcessor.BaseCommands(); // todo - make these cleaner, or by default...
	rpnProcessor.CommandsDone();

	

	if (argc <= 1)
	{
		ShowUsage(rpnProcessor);
		return -1;
	}
	bool verbose = false;

	State s;

	int argpos = 1; // skip initial exe
	while (argpos < argc && argv[argpos][0] == '-')
	{ // parse options
		string opt(argv[argpos]);
		argpos++;
		if (opt == "-s")
		{
			cout << "Executing script " << argv[argpos] << endl;
			rpnProcessor.GetScriptTokens(s.tokens, argv[argpos++]);
		}
		else if (opt == "--verbose")
		{
			cout << "Verbose = true\n";
			verbose = true;
		}
		else if (opt == "-a")
		{
			string t(argv[argpos++]);
			s.args.push_back(RPNLanguage::ToItem(t));
		}
		else {
			cerr << fmt::format("Unknown option {}\n", opt);
			return -1;
		}
	}
	// tokenize any other command line options
	for (int i = argpos; i < argc; ++i)
		s.tokens.emplace_back(argv[i]);

	cout << "Current path is " << fs::current_path() << '\n';
	const auto retval = rpnProcessor.Process(s, verbose) ? 1 : 0;
	cout << "Done\n";
	return retval;
}

