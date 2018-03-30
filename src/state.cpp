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
 * \brief State of the system
 *
 * Implementation of the methods of the class State to simulate
 * interacting active Brownian particles in dimension 2.
*/

#include <cmath>
#include <chrono>
// #include <iostream>
#include "state.h"

/*!
 * \brief Constructor of State
 *
 * Initializes the state of the system: particles randomly placed in a 2d box.
 *
 * \param _len Length of the box
 * \param _n_parts Number of particles
 * \param _pot_strength Strength of the interparticle potential
 * \param _temperature Temperature
 * \param _rot_dif Rotational diffusivity
 * \param _activity Activity
 * \param _dt Timestep
 */
State::State(const double _len, const long _n_parts,
	         const double _pot_strength, const double _temperature,
			 const double _rot_dif, const double _activity, const double _dt) :
	len(_len), n_parts(_n_parts), pot_strength(_pot_strength),
	activity(_activity), dt(_dt),
	// We seed the RNG with the current time
	rng(std::chrono::system_clock::now().time_since_epoch().count()),
	// Gaussian noise from the temperature
	noiseTemp(0.0, std::sqrt(2.0 * _temperature * dt)),
	// Gaussian noise from the rotational diffusivity
	noiseAngle(0.0, std::sqrt(2.0 * _rot_dif * dt)),
	boxes(_len, _n_parts)
{
	positions.resize(n_parts);
	angles.resize(n_parts);
	forces.resize(n_parts);

    std::uniform_real_distribution<double> rndPos(0, len);
    std::uniform_real_distribution<double> rndAngle(0, 2.0 * M_PI);

	for (long i = 0 ; i < n_parts ; ++i) {
		positions[i][0] = rndPos(rng);
		positions[i][1] = rndPos(rng);
		angles[i] = rndAngle(rng);
		forces[i][0] = 0;
		forces[i][1] = 0;
	}

	/* std::cout << boxes.getNBoxes() << std::endl;
	for (long i = 0 ; i < boxes.getNBoxes() ; ++i) {
		std::cout << "[" << i << "] :";
		for (auto it = boxes.getNbrsBegin(i) ; it != boxes.getNbrsEnd(i) ;
		     ++it) {
			std::cout << " " << *it;
		}
		std::cout << std::endl;
	} */
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation.
 */
void State::evolve() {
	double c, s;
	calcInternalForces();

	for (long i = 0 ; i < n_parts ; ++i) {
		// Computation of sin and cos
#ifdef __GNUC__
		sincos(angles[i], &s, &c);
#else
		s = std::sin(angles[i]);
		c = std::cos(angles[i]);
#endif
		// Internal forces +  Activity + Gaussian noise
		positions[i][0] += dt * (forces[i][0] + activity * c);
		positions[i][1] += dt * (forces[i][1] + activity * s);
		positions[i][0] += noiseTemp(rng);
		positions[i][1] += noiseTemp(rng); 
		// Rotational diffusion
		angles[i] += noiseAngle(rng);
	}

	enforcePBC();
}

/* \brief Compute the forces between the particles.
 *
 * Implement harmonic spheres.
 */
void State::calcInternalForces() {
    for (long i = 0 ; i < n_parts ; ++i) {
		forces[i][0] = 0;
		forces[i][1] = 0;
    }

	// Recompute the boxes
	boxes.update(&positions);

    for (long i = 0 ; i < n_parts ; ++i) {
		long k = boxes.getBoxOfPart(i);
		/*std::cout << "(" << positions[i][0] << ", " << positions[i][1]
			      << ") -> " << k << "\n";*/

		const auto it_b_b = boxes.getNbrsBegin(k);
		const auto it_b_e = boxes.getNbrsEnd(k);
		for(auto it_b = it_b_b ; it_b != it_b_e ; ++it_b) {
			const auto it_j_b = boxes.getPartsOfBoxBegin(*it_b);
			const auto it_j_e = boxes.getPartsOfBoxBegin(*it_b);
            for (auto it_j = it_j_b ; it_j != it_j_e && (*it_j) > i ; ++it_j) {
				double dx = positions[i][0] - positions[*it_j][0];
				double dy = positions[i][1] - positions[*it_j][1];
				// We want the periodized interval to be centered in 0
				pbcSym(dx, len);
				pbcSym(dy, len);
				double dr2 = dx * dx + dy * dy;

				if(dr2 < 1. && dr2 > 0.) {
					double u = pot_strength * (1.0 / std::sqrt(dr2) - 1.0);
					double fx = u * dx;
					double fy = u * dy;

					forces[i][0] += fx;
					forces[*it_j][0] -= fx;
					forces[i][1] += fy;
					forces[*it_j][1] -= fy;
				}
			}
		}
    }
}

/*
void State::calcInternalForces() {
    for (long i = 0 ; i < n_parts ; ++i) {
		forces[i][0] = 0;
		forces[i][1] = 0;
    }

    for (long i = 0 ; i < n_parts ; ++i) {
        for (long j = i + 1 ; j < n_parts ; ++j) {
			double dx = positions[i][0] - positions[j][0];
			double dy = positions[i][1] - positions[j][1];
			// We want the periodized interval to be centered in 0
			pbcSym(dx, len);
			pbcSym(dy, len);
			double dr2 = dx * dx + dy * dy;

            if(dr2 < 1. && dr2 > 0.) {
				double u = pot_strength * (1.0 / std::sqrt(dr2) - 1.0);
				double fx = u * dx;
				double fy = u * dy;

				forces[i][0] += fx;
				forces[j][0] -= fx;
				forces[i][1] += fy;
				forces[j][1] -= fy;
			}
        }
    }
}
*/

/* 
 * \brief Enforce periodic boundary conditions
 */
void State::enforcePBC() {
	for (long i = 0 ; i < n_parts ; ++i) {
		pbc(positions[i][0], len);
		pbc(positions[i][1], len);
		pbc(angles[i], 2.0 * M_PI);
	}
}
