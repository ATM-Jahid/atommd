#include <iostream>
#include <string>

void printSummary(std::string);

int main(int argc, char **argv) {
	std::string dot_in(argv[1]);

	printSummary(dot_in);
	std::cout << dot_in << '\n';
	return 0;
}

void printSummary(std::string dot_in) {
	std::string dot_out = dot_in.erase(dot_in.length()-2).append("out");
	std::cout << dot_out << '\n';
}
