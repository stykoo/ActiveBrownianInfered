/*
Copyright (C) Sorbonne Université (2018)
Contributor: Alexis Poncet <aponcet@lptmc.jussieu.fr>

This software (ActiveBrownian) intends to investigate the behavior
of interacting active Brownian particles in dimension 2.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \mainpage ActiveBrownian
 *
 * ActiveBrownian is a software to investigate the behavior of interacting
 * active Brownian particles (ABP) in dimension 2.
 * ABP are particles undergoing Brownian motion and having a motor
 * propulsing in a random direction with constant norm.
 */
/*!
 * \file main.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Main file of ActiveBrownian
*/

#include "simul.h"

/*!
 * \brief Main function
 *
 * Create and run the simulation.
 */
int main(int argc, char **argv) {
	Simul simulation(argc, argv);

	if (simulation.getStatus() == SIMUL_INIT_HELP) {
		return 0;
	} else if (simulation.getStatus() == SIMUL_INIT_FAILED) {
		return 1;
	}

	simulation.print();
	simulation.run();

	return 0;
}
