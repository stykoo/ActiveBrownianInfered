/*
Copyright (C) Sorbonne Université (2018)
Contributor: Alexis Poncet <aponcet@lptmc.jussieu.fr>

This file is part of ActiveBrownian.

ActiveBrownian is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ActiveBrownian is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ActiveBrownian.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file simul.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Simulation of the system
 *
 * Implementation of the methods of the class Simul to simulate
 * interacting active Brownian particles in dimension 2.
*/

#include <exception>
#include <thread>
#include <boost/program_options.hpp>
#include "simul.h"
#include "state.h"
#include "state3d.h"
#include "visu.h"
#include "visu3d.h"

namespace po = boost::program_options;

/*!
 * \brief Constructor of Simul
 *
 * Initializes the parameters of structure Simul
 * from the command-line arguments using boost::program_options.
 *
 * \param argc Number of arguments
 * \param argv Arguments
 */
Simul::Simul(int argc, char **argv) {
	status = SIMUL_INIT_SUCCESS;

	po::options_description opts("Options");
	opts.add_options()
		("rho,r", po::value<double>(&rho)->required(), "Density")
		("parts,n", po::value<long>(&n_parts)->required(),
		 "Number of particles")
		("eps,e", po::value<double>(&pot_strength)->required(),
		 "Strength of interparticle potential")
		("temp,T", po::value<double>(&temperature)->required(), "Temperature")
		("rdif,D", po::value<double>(&rot_dif)->required(),
		 "Rotational diffusivity")
		("activ,U", po::value<double>(&activity)->required(), "Activity")
		("dt,t", po::value<double>(&dt)->required(), "Timestep")
		("iters,I", po::value<long>(&n_iters)->required(),
		 "Number of time iterations")
		("3d,3", po::bool_switch(&sim3d), "Simulation in 3d instead of 2d")
		("sleep", po::value<int>(&sleep)->default_value(0),
		 "Number of milliseconds to sleep for between iterations")
		("help,h", "Print help message and exit")
		;

	try {
		po::variables_map vars;
		po::store(po::parse_command_line(argc, argv, opts), vars);

		// Display help and exit
		if (vars.count("help")) {
			std::cout << "Usage: " << argv[0] << " options\n";
			std::cout << opts << std::endl;
			status = SIMUL_INIT_HELP;
			return;
		}

        po::notify(vars);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		status = SIMUL_INIT_FAILED;
		return;
	}

	// Check if the values of the parameters are allowed
	if (notStrPositive(rho, "rho") || notStrPositive(n_parts, "n_parts")
		|| notPositive(pot_strength, "eps")
		|| notPositive(temperature, "T") || notPositive(rot_dif, "rot_dif")
		|| notPositive(activity, "actity") || notStrPositive(dt, "dt")
		|| notPositive(n_iters, "n_iters")) {
		status = SIMUL_INIT_FAILED;
		return;
	}

	if (sim3d) {
		len = std::cbrt(n_parts / rho);
	} else {
		len = std::sqrt(n_parts / rho);
	}
}

/*!
 * \brief Run the simulation
 *
 * Construct the state of the system and update it for the number
 * of iterations wanted. Also take care of launching the thread for
 * visualization.
 */
void Simul::run() {
	if (status != SIMUL_INIT_SUCCESS) {
		std::cerr << "You should not be runing a failed simulation..."
		          << std::endl;
		return;
	}

	if (sim3d) {
		// Initialize the state of the system
		State3d state(len, n_parts, pot_strength, temperature, rot_dif,
				      activity, dt);
		
		// Start thread for visualization
		Visu3d visu(&state, len, n_parts);
		std::thread thVisu(&Visu3d::run, &visu); 

		// Time evolution
		for (long t = 0 ; t < n_iters ; ++t) {
			state.evolve();
			if (sleep > 0) {
				std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
			}
		}

		thVisu.join();
	} else {
		// Initialize the state of the system
		State state(len, n_parts, pot_strength, temperature, rot_dif, activity,
					dt);
		
		// Start thread for visualization
		Visu visu(&state, len, n_parts);
		std::thread thVisu(&Visu::run, &visu); 

		// Time evolution
		for (long t = 0 ; t < n_iters ; ++t) {
			state.evolve();
			if (sleep > 0) {
				std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
			}
		}

		thVisu.join();
	}
}

/*!
 * \brief Print the parameters of the simulation
 */
void Simul::print() const {
	if (sim3d) {
		std::cout << "# [3d] ";
	} else {
		std::cout << "# [2d] ";
	}
	std::cout << "rho=" << rho << ", n_parts=" << n_parts
	          << ", pot_strength=" << pot_strength << ", temperature="
			  << temperature << ", rot_dif=" << rot_dif << ", activity="
			  << activity << ", dt=" << dt << ", nb_iters=" << n_iters << "\n";
	std::cout << std::endl;
}