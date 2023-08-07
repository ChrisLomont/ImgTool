#pragma once
// CSV file stuff
#include <memory>
#include <iostream>
#include <fstream>
#include "State.h"

	//{ "csvstart", " csvname header1 header2 ... headern n -> , start a CSV file with given headers", CsvOp},
	//{ "csvput"  , " val header csvname -> , stores val under header name in named csv",CsvOp },
	//{ "csvwrite", " csvname filename -> , ",CsvOp },
using namespace std;

namespace {
	class CSV
	{
	public:
		string name;
		vector<string> headers;
		unordered_map<string, vector<Item>> items;
	};
	using CSVPtr = shared_ptr<CSV>;

	unordered_map<string, CSVPtr> csvs;

	string ReadText(State & s)
	{ // do this way to allow numerical headers, filenames, etc
		const Item header = s.Pop();
		return FormatItem(header, false);
	}

}

inline void CsvOp(State& s, const string args)
	{
		if (args == "csvstart")
		{
			int n = s.PopInt();
			vector<string> headers;
			for (int i = 0; i < n; ++i)
			{
				string header = ReadText(s);
				headers.push_back(header);
			}
			auto csv = make_shared<CSV>();
			auto name = ReadText(s);
			csv->name = name;
			for (int i = n - 1; i >= 0; --i)
			{
				auto header = headers[i];
				csv->headers.push_back(header);
				csv->items[header] = {};
			}
			csvs[name] = csv;
		}
		else if (args == "csvput")
		{
			//{ "csvwrite", " csvname filename -> , ",CsvOp },
			auto header = ReadText(s);
			auto csvname = ReadText(s);
			Item item = s.Pop();
			if (!contains(csvs, csvname))
			{
				throw runtime_error(fmt::format("No CSV named {}", csvname));
			}
			auto& items = csvs[csvname]->items;
			if (contains(items, header))
				items[header].push_back(item);
			else
				throw runtime_error(fmt::format("CSV {} does not have header {}", csvname, header));
		}
		else if (args == "csvwrite")
		{
			//{ "csvwrite", " csvname filename -> , ",CsvOp },
			auto filename = ReadText(s);
			auto csvname = ReadText(s);
			const auto & csv = csvs[csvname];
			
			ofstream outfile;
			outfile.open(filename);
			// write headers
			for (auto h : csv->headers)
				outfile << h << ", ";
			outfile << endl;
			
			// find max row count
			int maxrow = 0;
			for (auto h : csv->headers)
			{
				const auto& v = csv->items[h];
				const int size = static_cast<int>(v.size());
				maxrow = max(maxrow, size);
			}

			// write rows
			int row = 0;
			while (row < maxrow)
			{
				for (auto h : csv->headers)
				{
					const auto& v = csv->items[h];					
					if (v.size() <= row)
					{
						outfile << "<null>, ";
					}
					else
					{
						outfile << FormatItem(v[row],false) << ", ";
					}
				}
				outfile << endl;
				++row;
			}

			outfile.close();
			}	
		else
			throw runtime_error("");
	}
