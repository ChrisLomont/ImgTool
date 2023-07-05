#pragma once

#include <chrono>

using namespace std; // todo remove

class Timer
{
public:

	using clock_t = std::chrono::high_resolution_clock;
	using duration_t = clock_t::duration;

private:
	clock_t::time_point start_time;
	clock_t::time_point split_time;
public:

	Timer() { reset(); }

	// total elapsed time
	duration_t get_elapsed_time() const { return (clock_t::now() - start_time); }
	// call to get this split time, starts new split
	duration_t get_split_time() {
		const auto c = clock_t::now();
		const auto e = c - split_time;
		split_time = c;
		return e;
	}
	// set elapsed and split to 0
	void reset() { split_time = start_time = clock_t::now(); }

	// use as cout << Timer::format_us(duration) << endl;
	static string format_us(const duration_t& duration) {
		using namespace std::chrono;
		double t = duration_cast<microseconds>(duration).count();
		return format("{:0.4f}us", t);
	}
};