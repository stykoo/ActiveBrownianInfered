/*
Copyright (C) Sorbonne Universit√© (2018)
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
#include <boost/program_options.hpp>
#include "observables.h"
#include "simul.h"
#include "state.h"

#ifndef NOVISU
#include <thread>
#include "visu/visu.h"
#endif

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
		("temp,T", po::value<double>(&temperature)->required(), "Temperature")
		("rdif,D", po::value<double>(&rot_dif)->required(),
		 "Rotational diffusivity")
		("activ,U", po::value<double>(&activity)->required(), "Activity")
		("dt,t", po::value<double>(&dt)->required(), "Timestep")
		("iters,I", po::value<long>(&n_iters)->required(),
		 "Number of time iterations")
		("itersTh,J", po::value<long>(&n_iters_th)->default_value(0),
		 "Number of time iterations of thermalization")
		("skip,S", po::value<long>(&skip)->default_value(100),
		 "Iterations between two computations of observables")
		("nfuns", po::value<long>(&n_funs)->required(),
		 "Number of radial functions")
		("nmodes", po::value<long>(&n_modes)->required(),
		 "Number of angular modes")
		("infered",
		 po::value<std::string>(&fname_infered)->default_value("infered"),
		 "Filename for infered coefficients")
		("output,O",
		 po::value<std::string>(&output)->default_value("observables.h5"),
		 "Name of the output file")
		("stepr,s",
		 po::value<double>(&step_r)->default_value(0.2),
		 "Spatial resolution for correlations")
		("divAngle,d",
		 po::value<long>(&n_div_angle)->default_value(40))
#ifndef NOVISU
		("sleep", po::value<int>(&sleep)->default_value(0),
		 "Number of milliseconds to sleep for between iterations")
#endif
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
		|| notPositive(temperature, "T") || notPositive(rot_dif, "rot_dif")
		|| notPositive(activity, "actity") || notStrPositive(dt, "dt")
		|| notPositive(n_iters, "n_iters")) {
		status = SIMUL_INIT_FAILED;
		return;
	}

	len = std::sqrt(n_parts / rho);
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
	
	// Initialize the state of the system
	Infered infered(n_funs, n_modes, fname_infered + "_r.dat",
			        fname_infered + "_t.dat", fname_infered + "_o.dat");
	State state(len, n_parts, temperature, rot_dif, activity,
				dt, infered);
	Observables obs(len, n_parts, step_r, n_div_angle);
	obs.computeForces(infered);
	
#ifndef NOVISU
	// Start thread for visualization
	Visu visu(&state, len, n_parts);
	std::thread thVisu(&Visu::run, &visu); 
#endif

	// Thermalization
	for (long t = 0 ; t < n_iters_th ; ++t) {
		state.evolve();
	}
	// Time evolution
	for (long t = 0 ; t < n_iters ; ++t) {
		state.evolve();
		if (t % skip == 0) {
			obs.compute(&state);
		}
#ifndef NOVISU
		if (sleep > 0) {
			std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
		}
#endif
	}

	obs.writeH5(output, rho, n_parts, temperature, rot_dif,
				activity, dt, n_iters, n_iters_th, skip);
	//state.dump();

#ifndef NOVISU
		thVisu.join();
#endif
}

/*!
 * \brief Print the parameters of the simulation
 */
void Simul::print() const {
	std::cout << "# [2d] ";
	std::cout << "rho=" << rho << ", n_parts=" << n_parts
	          << ", temperature="
			  << temperature << ", rot_dif=" << rot_dif << ", activity="
			  << activity << ", dt=" << dt << ", n_iters=" << n_iters
			  << ", n_iters_th=" << n_iters_th << ", skip=" << skip
			  << "\n";
	std::cout << std::endl;
}
