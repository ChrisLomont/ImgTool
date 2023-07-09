#pragma once
#include <string>
#include <algorithm>
#include <regex>

using namespace std; // todo - remove

/*----------------- Utils -------------------------------*/
//namespace {
	string ToUpper(const string& text)
	{
		string s(text);
		transform(s.begin(), s.end(), s.begin(), ::toupper);
		return s;
	}
	string ToLower(const string& text)
	{
		string s(text);
		transform(s.begin(), s.end(), s.begin(), [](char c) {return std::tolower(c); });
		return s;
	}
	double ParseDouble(const string& item)
	{
		// todo - error checking?
		return stod(item);
	}

	// can parse to double? (requires decimal?)
	bool IsDouble(const string& text)
	{
		return regex_match(text, regex("^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$"));		
	}

	const std::string WHITESPACE = " \n\r\t\f\v";

	std::string ltrim(const std::string& s)
	{
		size_t start = s.find_first_not_of(WHITESPACE);
		return (start == std::string::npos) ? "" : s.substr(start);
	}

	std::string rtrim(const std::string& s)
	{
		size_t end = s.find_last_not_of(WHITESPACE);
		return (end == std::string::npos) ? "" : s.substr(0, end + 1);
	}

	std::string trim(const std::string& s) {
		return rtrim(ltrim(s));
	}

// }