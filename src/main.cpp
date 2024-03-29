#include <cassert>
#include "simul.h"

int main(int argc, char **argv) {
	assert(argc > 2);
	Simul simulation(argv[1]);

	if (simulation.getStatus() == SIMUL_INIT_FAILED)
		return 1;

	simulation.print();
	simulation.run();
	simulation.save(argv[2]);

	return 0;
}
