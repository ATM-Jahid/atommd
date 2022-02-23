#include <iostream>
#include <stdlib.h>
using namespace std;

int main() {
	srand(17);
	cout << rand() % 1000 << '\n';
	cout << (rand() % 1000)/999.9 << '\n';

	return 0;
}
