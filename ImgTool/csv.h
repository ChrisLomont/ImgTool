#pragma once
// CSV file stuff
#include <memory>
#include "State.h"

	//{"csvstart", " csvname header1 header2 ... n -> , start a CSV file with given headers", CsvOp},
	//{ "csvput"  , " val header csvname -> , stores val under header name in named csv",CsvOp },
	//{ "csvwrite", " csvname filename -> , ",CsvOp },

class CSV
{
public:
	vector<string> headers;
	unordered_map<string, vector<Item>> items;
};
using CSVPtr = shared_ptr<CSV>;

unordered_map<string, CSVPtr> csvs;

void CsvOp(State& s, const string args)
	{
		if (args == "csvstart")
		{
			int n = s.PopInt();
			vector<string> headers;
			for (int i = 0; i < n; ++i)
			{
				string header = s.Pop<string>();
				headers.push_back(header);
			}
			string csvname = s.Pop<string>();
		}
		else if (args == "csvput")
		{

		}
		else if (args == "csvwrite")
		{

		}
		throw runtime_error("");
	}
