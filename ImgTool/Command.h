#pragma once

#include <functional>
#include <string>
#include "State.h"

using namespace std;

using Op = function<void(State&, const string& args)>;

struct Command {
	string name;
	string description;
	Op op;
};