#include <iostream>
#include <random>

int main() {
	int Nexp = 10000;
	int Star = 100;

	std::default_random_engine rand_gen;
	std::normal_distribution<double> normal_dist(10.0, 3.0);

	int p[20];

	for (int i = 0; i < Nexp; i++) {
		double number = normal_dist(rand_gen);
		if (number >= 0 && number < 20) {
			p[int(number)]++;
		}
	}

	for (int i = 0; i < 20; i++) {
		std::cout << i << '-' << i+1 << ": ";
		std::cout << std::string(p[i]*Star/Nexp, '*') << '\n';
	}

	return 0;
}
