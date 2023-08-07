#pragma once
#include <unordered_map>
#include <stack>
#include "Stack.h"

using namespace std; // todo - remove

class State : public Stack
{
	unordered_map<string, int> labels;
	unordered_map<string, Item> vars;

	struct Loop
	{
		int min, max, loopPosition;
		vector<Item> items;
		int index{ 0 }; // number sent out
	};
	vector<Loop> loops;
	int nextLoopIndex = 0;
	void AddLoop(int a, int b, vector<Item> & items)
	{
		while (loops.size() <= nextLoopIndex)
			loops.push_back(Loop{});
		loops[nextLoopIndex++] = Loop{a,b,programPosition,items};
	}
	void DoLoop()
	{
		if (nextLoopIndex <= 0)
			throw runtime_error("No more loops to do");
		auto& loop = loops[nextLoopIndex - 1];

		// how many to do total
		const int numItems = abs(loop.min - loop.max) + 1;

		// see if done, if not, push items, else remove loop
		if (loop.index < numItems)
		{
			// if items, push one (backwards, makes easier to use)
			if (loop.index < loop.items.size())
				Push(loop.items[loop.items.size()-1-loop.index]);

			// push index
			const int di = loop.max > loop.min ? 1 : -1;
			const int index = loop.min + loop.index * (di);
			Push(index);

			loop.index++;

			// goto loop top
			programPosition = loop.loopPosition;
		}
		else
		{ // remove loop
			nextLoopIndex--;
		}
	}

	std::stack<int> programCounterStack;	

	Timer timer{};

public:
	vector<string> tokens;

	State() {
		randState = randSeed(42);
	}

	vector<Item> args;

	// elapsed us time
	int64_t elapsedUs()
	{
		const auto duration = timer.get_elapsed_time();
		return std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
	}



	//  0 = none, 1 = info, 2 = all
	int verbosity = 1;
	int exitCode{ 0 }; // exit code on halt
	bool inSubroutineDefinition{ false };

	uint64_t randState{ 1234 };

	void StateOp(const string & args)
	{
		if (args == "label")
		{
			auto label = Pop<string>();
			labels[label] = programPosition;
		}
		else if (args == "subroutine")
		{
			auto label = Pop<string>();
			labels[label] = programPosition;
			inSubroutineDefinition = true;
		}
		else if (args == "endsub")
		{
			inSubroutineDefinition = false;
		}
		else if (args == "gosub")
		{
			auto label = Pop<string>();
			programCounterStack.push(programPosition);
			programPosition = labels[label];
		}
		else if (args == "return")
		{
			programPosition = programCounterStack.top();
			programCounterStack.pop();
		}
		else if (args == "ja")
		{
			auto label = Pop<string>();
			programPosition = labels[label];
		}
		else if (
			args == "je" ||
			args == "jne" ||
			args == "jlt" ||
			args == "jle" ||
			args == "jgt" ||
			args == "jge"
			)
		{
			auto label = Pop<string>();
			auto item2 = Pop<double>();
			auto item1 = Pop<double>();

			bool test1 = args == "je" && item1 == item2;
			bool test2 = args == "jne" && item1 != item2;
			bool test3 = args == "jlt" && item1 < item2;
			bool test4 = args == "jle" && item1 <= item2;
			bool test5 = args == "jgt" && item1 > item2;
			bool test6 = args == "jge" && item1 >= item2;

			if (test1 || test2 || test3 || test4 || test5 || test6)
			{
				programPosition = labels[label];
			}
		}
		else if (args == "halt")
		{
			exitCode = PopInt();
			programPosition = 1 << 30; // out of bounds - todo - make const?, flag?
		}
		else if (args == "sto")
		{
			auto name = Pop<string>();
			auto item = Pop();
			vars[name] = item;
		}
		else if (args == "rcl")
		{
			auto name = Pop<string>();
			if (!vars.contains(name))
			{
				throw runtime_error(fmt::format("Cannot 'rcl' variable {}", name));
			}
			auto item = vars[name];
			Push(item);
		}
		else if (args == "dumpstate")
		{
			cout << "State vars:\n";
			for (auto& i : vars)
			{
				cout << fmt::format("{}: {}\n", i.first, FormatItem(i.second, true));
			}
			cout << "State labels:\n";
			for (auto& i : labels)
			{
				cout << fmt::format("{}: {}\n", i.first, i.second);

			}
		}
		else if (args == "rangeloop")
		{
			auto b = PopInt();
			auto a = PopInt();
			vector<Item> t;
			AddLoop(a,b,t);
			DoLoop();
		}
		else if (args == "itemloop")
		{
			auto n = PopInt();
			vector<Item> items;
			for (auto i = 0; i < n; ++i)
				items.push_back(Pop());
			AddLoop(0,n-1,items);
			DoLoop();
		}
		else if (args == "endloop")
		{
			DoLoop();
		}
		else if (args == "system")
		{
			auto cmd = Pop<string>();
			auto retval = system(cmd.c_str());
			Push(retval);
		}
		else if (args=="verbosity")
		{
			auto v = PopInt();
			verbosity = v;
		}
		else if (args == "if")
		{
			auto isTrue = Pop<double>() != 0;

			auto m = PopInt();
			auto n = PopInt();
			vector<Item> f = PopN(m);
			vector<Item> t = PopN(n);
			if (isTrue)
				PushN(t);
			else
				PushN(f);
		}
		else if (args == "->str")
		{
			auto v = Pop();
			Push(FormatItem(v, false));
		}
		else 
			throw runtime_error(fmt::format("unknown state op {}",args));
	}
	int programPosition{0};

};