#pragma once
#include <unordered_map>
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
		int numItems = abs(loop.min - loop.max) + 1;

		// see if done, if not, push items, else remove loop
		if (loop.index < numItems)
		{
			// if items, push one (backwards, makes easier to use)
			if (loop.index < loop.items.size())
				Push(loop.items[loop.items.size()-1-loop.index]);

			// push index
			int di = loop.max > loop.min ? 1 : -1;
			int index = loop.min + loop.index * (di);
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

public:
	//  0 = none, 1 = info, 2 = all
	int verbosity = 1;

	void StateOp(const string & args)
	{
		if (args == "label")
		{
			auto label = Pop<string>();
			labels[label] = programPosition;
		}
		else if (args == "ja")
		{
			auto label = Pop<string>();
			programPosition = labels[label];
		}
		else if (args == "je")
		{
			auto label = Pop<string>();
			auto v = Pop<double>();
			if (v != 0)
			{
				programPosition = labels[label];
			}
		}
		else if (args == "halt")
		{
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
			auto item = vars[name];
			Push(item);
		}
		else if (args == "dumpstate")
		{
			cout << "State vars:\n";
			for (auto& i : vars)
			{
				cout << format("{}: {}\n", i.first, FormatItem(i.second, true));
			}
			cout << "State labels:\n";
			for (auto& i : labels)
			{
				cout << format("{}: {}\n", i.first, i.second);

			}
		}
		else if (args == "rangeloop")
		{
			auto b = ParseInt(Pop<string>());
			auto a = ParseInt(Pop<string>());
			vector<Item> t;
			AddLoop(a,b,t);
			DoLoop();
		}
		else if (args == "itemloop")
		{
			auto n = ParseInt(Pop<string>());
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
			Push(format("{}",retval));
		}
		else if (args=="verbosity")
		{
			auto v = ParseInt(Pop<string>());
			verbosity = v;
		}
		else if (args == "if")
		{
			auto v = ParseDouble(Pop<string>());
			auto m = ParseInt(Pop<string>());
			auto n = ParseInt(Pop<string>());
			vector<Item> f = PopN(m);
			vector<Item> t = PopN(n);
			if (v != 0)
				PushN(t);
			else
				PushN(f);
		}
		else 
			throw runtime_error(format("unknown state op {}",args));
	}
	int programPosition{0};

};