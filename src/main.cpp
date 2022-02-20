#include <iostream>

#include "types.hpp"
#include "vec_cal.hpp"

int main(int argc, char **argv) {
	std::cout << "hello world\n";

	vecR u, v, mul, div;

	vecSet(u, 10, 20);
	vecSet(v, -3, -9);
	vecMul(mul, u, v);
	vecDiv(div, u, v);

	std::cout << mul.x << div.y << '\n';
	vecSelfScaleAdd(u, 3, v);
	std::cout << u.x << u.y << '\n';

	vecR p, q;
	vecSetAll(p, -12);
	vecSetZero(q);

	vecR reg;
	vecSetAll(reg, 20);
	vecWrapAll(p, reg);

	std::cout << p.x << p.y << '\n';

	return 0;
}
