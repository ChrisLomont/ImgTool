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
			exitCode = Pos() > 0 ? PopInt() : 1; // default
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
			auto m = PopInt();
			auto n = PopInt();
			vector<Item> f = PopN(m);
			vector<Item> t = PopN(n);
			auto isTrue = Pop<double>() != 0;
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
		else if (args == "->list")
		{
			// { "->list"," item1 item2 ... itemn n  -> list of items, convert n items into a list", StateOp },
		
			auto n = PopInt();
			Item it;
			for (auto i = 0; i < n; ++i)
				it.arr.push_back(Pop());
			it.type = ItemType::List;
			Push(it);
		}
		else if (args == "list->")
		{
			//{ "list->"," list -> item1 item2 ... itemn n, list of n items out", StateOp },
			auto lst = PopList();
			auto n = (int)lst.size();
			for (auto& it : lst)
				Push(it);
			Push(n);

		}
		else if (args == "listlen")
		{
			//{ "listlen"," list -> list list_length, get length of list", StateOp },
			auto it = Pop();
			if (!holds_alternative<std::vector<Item>>(it))
				throw runtime_error("Expected list");
			auto n = (int)it.arr.size();
			Push(it);
			Push(n);
		}
		else if (args == "listget")
		{
			// { "listget", " list k -> list item_k, get kth item from list, 0 indexed", StateOp },
			throw runtime_error("listget not implemented");

		}
		else if (args == "listset")
		{
			// { "listset"," list item k -> list, set kth item from list, 0 indexed", StateOp },
			throw runtime_error("listset not implemented");
		}
		else if (args == "sublist")
		{
			// { "sublist"," list a b -> sublist, get sublist of items a (inclusive) to b (exclusive) 0 indexed", StateOp },
			throw runtime_error("sublist not implemented");
		}
		else if (args == "listins")
		{
			// { "listins"," list item k -> list, insert item at index k, 0 indexed", StateOp },
			throw runtime_error("listins not implemented");
			}
		else if (args == "listdel")
		{
			// { "listdel"," list k -> list, delete the k item, 0 indexed", StateOp },
			throw runtime_error("listdel not implemented");
		}
		else if (args == "listappend")
		{
			// { "listappend"," list item -> append item to list", StateOp },
			throw runtime_error("listappend not implemented");
		}
		else if (args == "listjoin")
		{
			// { "listjoin"," list1 list2 -> list, join lists 1 and 2", StateOp },
			throw runtime_error("listjoin not implemented");
		}
		else 
			throw runtime_error(fmt::format("unknown state op {}",args));
	}
	int programPosition{0};

};