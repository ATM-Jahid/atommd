#include <iostream>
#define CELL_LIST
#define NEIGH_LIST

int main() {
#ifdef CELL_LIST
	std::cout << "hello cells!\n";
#endif

#ifdef NEIGH_LIST
	std::cout << "hello neighbors!\n";
#endif

#if defined(CELL_LIST) && defined(NEIGH_LIST)
	std::cout << "hello cell neighbors!\n";
#endif

	return 0;
}
