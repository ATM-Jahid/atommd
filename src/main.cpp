#include <iostream>

#include "types.hpp"
#include "vec_cal.hpp"

int main(int argc, char **argv) {
	std::cout << "hello world\n";

	vecR z;
	vecSetAll(z, 1);

	std::cout << z.x << '\n';

	return 0;
}
