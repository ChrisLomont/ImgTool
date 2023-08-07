#include <exception>
#include <string>
// exception to use in system
class image_exception : public std::exception
{
public:
	std::string msg;
	image_exception(const std::string& msg) : msg(msg) {}

	const char* what() {
		return msg.c_str();
	}
};