#pragma once
#include <string>
#include <algorithm>
#include <regex>
#include <unordered_map>

using namespace std; // todo - remove

// math utility:

// sign of value -1,0,1
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
} // https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c


/*----------------- Utils -------------------------------*/
//namespace {
inline string ToUpper(const string& text)
	{
		string s(text);
		transform(s.begin(), s.end(), s.begin(), ::toupper);
		return s;
	}
inline string ToLower(const string& text)
	{
		string s(text);
		transform(s.begin(), s.end(), s.begin(), [](char c) {return std::tolower(c); });
		return s;
	}
inline double ParseDouble(const string& item)
	{
		// todo - error checking?
		return stod(item);
	}

	// can parse to double? (requires decimal?)
inline bool IsDouble(const string& text)
	{
		return regex_match(text, regex("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$"));		
	}

	const std::string WHITESPACE = " \n\r\t\f\v";

	inline std::string ltrim(const std::string& s)
	{
		const size_t start = s.find_first_not_of(WHITESPACE);
		return (start == std::string::npos) ? "" : s.substr(start);
	}

	inline std::string rtrim(const std::string& s)
	{
		const size_t end = s.find_last_not_of(WHITESPACE);
		return (end == std::string::npos) ? "" : s.substr(0, end + 1);
	}

	inline std::string trim(const std::string& s) {
		return rtrim(ltrim(s));
	}

	template<typename TKey, typename TValue>
	bool contains(const std::unordered_map<TKey,TValue> & map, const TKey & key)
	{
		return map.find(key) != map.end();
	}

	

// }