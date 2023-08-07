#pragma once
#include <cstdint>



/*
with default seed set to 42, these are first few 64 bit values
0x27a53829edf003a9,
0xdf28458e5c04c31c,
0x2756dc550bc36037,
0xa10325553eb09ee9,
0x40a0fccb8d9df09f,
0x5c2047cfefb5e9ca
*/

// Random 64 bit unsigned integer, uniformly distributed, period 2^64
// algorithm PCG RXS M XS 64 (LCG) 
// has period 2^64, state 64 bits, speed 1.01 ns/rand (from paper setup)
uint64_t rand64(uint64_t& state)
{
	const uint64_t PcgDefaultMultiplier64 = 6364136223846793005ULL;
	const uint64_t PcgDefaultIncrement64 = 1442695040888963407ULL;
	const auto temp = ((state >> (static_cast<int>(state >> 59) + 5)) ^ state) * 12605985483714917081ULL;
	state = state * PcgDefaultMultiplier64 + PcgDefaultIncrement64;
	return (temp >> 43) ^ temp;
}

// set rand seed, return state
// should be called to get a state
uint64_t randSeed(uint64_t seed = 42)
{ // set seed according to paper
	uint64_t state = 0;
	rand64(state);
	state += seed;
	rand64(state);
	return state;
}

// uniformly distributed 32 bit random
uint32_t rand32(uint64_t& state)
{
	return static_cast<uint32_t>(rand64(state) >> 20);
}


/// <summary>
 /// Uniform random integer in [min,max)
 /// </summary>
 /// <param type="max"></param>
 /// <returns></returns>
int32_t randUniform(uint64_t& state, int32_t min, int32_t max)
{
	if (max <= min) return min;

	const int delta = max - min;

	if (delta <= 0) return 0;
	const auto threshold = static_cast<uint32_t>((1ULL << 32) - (uint64_t)delta) % delta;

	while (true)
	{
		const auto r = rand32(state);
		if (r >= threshold)
			return static_cast<int32_t>(r % delta) + min;
	}
}
