#include <iostream>

#include "types.hpp"
#include "vec_cal.hpp"

int main(int argc, char **argv) {
	std::cout << "hello world\n";

	int a = 3, b = 4, sum, diff;
	//vecAdd(sum, a, b);
	std::cout << a << b << sum << '\n';
	vecAdd(sum, a, b);
	std::cout << a << b << sum << '\n';
	vecSub(diff, a, b);
	std::cout << a << b << diff << '\n';

	Prop abc;
	abc.val = 20;
	abc.sum = 40;
	abc.sum2 = -2;

	std::cout << abc.val << abc.sum2 << '\n';

	return 0;
}
