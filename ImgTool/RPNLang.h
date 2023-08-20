#pragma once
#include <string>
#include <regex>
#include <filesystem>
#include <set>
#include <fstream>
#include <vector>


// some compilers puke (clang!), even with c++ 20, so replace <format> with older one
//#include <format>
#include "fmt/fmt/format.h"


#include "State.h"
#include "Stack.h"
#include "Command.h"
#include "Utils.h"


using namespace std; // todo - refactor this out
namespace fs = std::filesystem;

// todo - refactor this, clean up

// RPN Language
class RPNLanguage
{

	// todo - make non-static
	static inline int VERSION_MAJOR{0};
	static inline int VERSION_MINOR{0};
public:

	// set version to report to code
	RPNLanguage(int majorVersion, int minorVersion)
	{
		VERSION_MAJOR = majorVersion;
		VERSION_MINOR = minorVersion;
	}


	// todo - clean up command list and command handling
	// where all commands are aggregated, ordered for usage messages
	vector<Command> commandList = {
	};


	// todo - make non-static
	static void StateOp(State& s, const string& args)
	{
		s.StateOp(args);
	}

	// todo - make non-static
	static void StackOp(State& s, const string& args)
	{
		if (s.verbosity >= 1)
			cout << fmt::format("Stack operation: {}\n", args);

		s.StackOp(args);
	}

	// todo - make this auto before any execution
	void CommandsDone()
	{
		for (auto& c : commandList)
		{
			commandMap[c.name] = c;
		}
	}


	// todo - make non-static
	static void RandOp(State& s, const string& args)
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
		else throw runtime_error(fmt::format("Unknown command {}", args));
	}

	// todo - make non-static
	static void Print(State& s, const string& args)
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
				cout << fmt::format("{}: ", i);
			cout << FormatItem(v, s.verbosity >= 2) << " ";
		}
		cout << endl;
	}

	static void GetFiles(State& s)
	{
		const auto regtext = s.Pop<string>();
		const auto path = s.Pop<string>();
		vector<string> files;
		const regex reg(regtext);

		for (const auto& entry : fs::directory_iterator(path))
		{
			if (entry.is_regular_file())
			{
				const auto& path1 = entry.path();
				const auto fn = path1.filename().string();
				if (regex_match(fn, reg))
					files.push_back(path1.string());
			}
		}
		sort(files.begin(), files.end());
		for (const auto& f : files)
		{
			s.Push(f);
		}
		s.Push(static_cast<int>(files.size()));
	}

	// todo - make this and other needed items non-static
	inline static set<string> filesRead;

	static void GetScriptTokens(vector<string>& tokens, const string& filename)
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
				//if (token[0] == '"') token = token.substr(1, token.size() - 2); // string keeps quotes, stripped in execution
				//cout << format(" => <{}>\n", token);
				tokens.push_back(token);
			}
		}
	}

	// todo - make non-static
	static void SystemOp(State& s, const string& args)
	{
		if (args == "getfiles")
		{
			GetFiles(s);
		}
		else if (args == "include")
		{
			//	{"include", " filename -> , include file as text, each file included at most once", IncludeFile},

			const auto filename = s.Pop<string>();

			vector<string> tokens;
			GetScriptTokens(tokens, filename);
			s.tokens.insert(next(s.tokens.begin(), s.programPosition), tokens.begin(), tokens.end());
		}
		else if (args == "version")
		{
			// version
			s.Push(VERSION_MAJOR);
			s.Push(VERSION_MINOR);
		}
		else if (args == "timeus")
		{
			s.Push(static_cast<double>(s.elapsedUs()));
		}
		else if (args == "arg")
		{
			const auto n = s.PopInt();
			s.Push(s.args[n]);
		}
		else if (args == "argcount")
		{
			s.Push(static_cast<int>(s.args.size()));
		}
		else throw runtime_error("Unknown system op");
	}

	// command map for quick lookup
	unordered_map<string, Command> commandMap;

	void BaseCommands()
	{
		vector<Command> systemCommands = {
			// system
			{"files","path regex -> f1 f2 ... fn n, reads files matching regex, pushes on stack with count",SystemOp},
			{"version"," -> major minor, get version", SystemOp},
			{"timeus"," -> time_us, get elapsed time in us", SystemOp},
			{"rand","a b -> rand32(a,b), uniform random integer in [a,b)",RandOp},
			{"srand","seed -> , set random seed to integer seed",RandOp},
			{"arg", " n -> arg, get command line arg n, passed via -a item, n = 1,2,...",SystemOp},
			{"argcount", "  -> argcount, count of command line args passed via via -a",SystemOp},

			{"include", " filename -> , include file as text, each file included at most once",SystemOp},
		};

		vector<Command> stackCommands = {
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
		};

		vector<Command> printCommands = {
			// printing
			{"print","item -> , prints top item",Print},
			{"printn","x1 x2 ... xn  n -> , prints top N items",Print},
			{"endl", "-> endline, pushes an endline string", Print},
		};

		vector<Command> stateCommands = {
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

			// list
			{ "->list"," item1 item2 ... itemn n  -> list of items, convert n items into a list", StateOp },
			{ "list->"," list -> item1 item2 ... itemn n, list of n items out", StateOp },
			{ "listlen"," list -> list list_length, get length of list", StateOp },
			{ "listget"," list k -> list item_k, get kth item from list, 0 indexed", StateOp },
			{ "listset"," list item k -> list, set kth item from list, 0 indexed", StateOp },
			{ "sublist"," list a b -> sublist, get sublist of items a (inclusive) to b (exclusive) 0 indexed", StateOp },
			{ "listins"," list item k -> list, insert item at index k, 0 indexed", StateOp },
			{ "listdel"," list k -> list, delete the k item, 0 indexed", StateOp },
			{ "listappend"," list item -> append item to list", StateOp },
			{ "listjoin"," list1 list2 -> list, join lists 1 and 2", StateOp },
		};

		commandList.insert(commandList.end(), systemCommands.begin(), systemCommands.end());
		commandList.insert(commandList.end(), stackCommands.begin(), stackCommands.end());
		commandList.insert(commandList.end(), printCommands.begin(), printCommands.end());
		commandList.insert(commandList.end(), stateCommands.begin(), stateCommands.end());

	}


	bool Process(State& state, bool verbose)
	{
		state.verbosity = verbose ? 2 : 1;

		state.programPosition = 0;
		try {
			while (state.programPosition < state.tokens.size())
			{
				const string& token = state.tokens[state.programPosition];
				state.programPosition++; // next position

				if (state.inSubroutineDefinition && token != "endsub")
					continue; // do nothing

				auto it = commandMap.find(token);
				if (it != commandMap.end())
				{ // execute operation
					const auto& c = it->second;
					if (state.verbosity >= 2)
						cout << fmt::format("Executing {}\n", c.name);
					c.op(state, token);
				}
				else {
					// was no item, push it
					// strip outer string quotes here 
					auto s = token; // get copy
					if (s[0] == '"') s = s.substr(1, s.size() - 2); // string keeps quotes, stripped in execution
					if (state.verbosity == 2)
						cout << fmt::format("Pushing {}\n", s);
					state.Push(ToItem(s));
				}

			}
		}
		catch (runtime_error& e)
		{
			cerr << fmt::format("Exception {}\n", e.what());
			return false;
		}
		catch (...)
		{
			cerr << fmt::format("Exception: should not reach this handler!");
			return false;
		}
		cout << fmt::format("Stack depth on exit {}, exit code {}\n", state.Pos(),state.exitCode);
		return true;
	}

	static Item ToItem(const string& text)
	{
		if (IsDouble(text))
			return Item(ParseDouble(text));
		return Item(text);
	}

	
};