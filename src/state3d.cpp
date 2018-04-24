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
 * \brief State of the system in dimension 3
 *
 * Implementation of the methods of the class State to simulate
 * interacting active Brownian particles in dimension 3.
*/

#include <cmath>
#include <chrono>
// #include <iostream>
#include "state.h"
#include "state3d.h"

/*!
 * \brief Constructor of State3d
 *
 * Initializes the state of the system: particles randomly placed in a 3d box.
 *
 * \param _len Length of the box
 * \param _n_parts Number of particles
 * \param _pot_strength Strength of the interparticle potential
 * \param _temperature Temperature
 * \param _rot_dif Rotational diffusivity
 * \param _activity Activity
 * \param _dt Timestep
 */
State3d::State3d(const double _len, const long _n_parts,
	             const double _pot_strength, const double _temperature,
			     const double _rot_dif, const double _activity,
				 const double _dt) :
	len(_len), n_parts(_n_parts), pot_strength(_pot_strength),
	activity(_activity), dt(_dt),
	// We seed the RNG with the current time
	rng(std::chrono::system_clock::now().time_since_epoch().count()),
	// Gaussian noise from the temperature
	noiseTemp(0.0, std::sqrt(2.0 * _temperature * dt)),
	// Standard deviation of gaussian noise from the rotational diffusivity
	stddevOrient(std::sqrt(2.0 * _rot_dif * dt)),
	boxes(_len, _n_parts)
{
	positions.resize(n_parts);
	// Automatic random intialization
	for (long i = 0 ; i < n_parts ; ++i) {
		orients.emplace_back(rng);
	}
	forces.resize(n_parts);

    std::uniform_real_distribution<double> rndPos(0, len);

	for (long i = 0 ; i < n_parts ; ++i) {
		positions[i][0] = rndPos(rng);
		positions[i][1] = rndPos(rng);
		positions[i][2] = rndPos(rng);
		forces[i][0] = 0;
		forces[i][1] = 0;
		forces[i][2] = 0;
	}

	/* std::cout << boxes.getNBoxes() << std::endl;
	for (long i = 0 ; i < boxes.getNBoxes() ; ++i) {
		std::cout << "[" << i << "] :";
		for (auto it = boxes.getNbrsBegin(i) ; it != boxes.getNbrsEnd(i) ;
		     ++it) {
			std::cout << " " << *it;
		}
		std::cout << std::endl;
	} 
	std::cout << _temperature << " " << stddevOrient << " " << activity
		<< std::endl; */
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation.
 */
void State3d::evolve() {
	calcInternalForces();

	for (long i = 0 ; i < n_parts ; ++i) {
		// Internal forces +  Activity + Gaussian noise
		positions[i][0] += dt * (forces[i][0] + activity * orients[i].getX());
		positions[i][1] += dt * (forces[i][1] + activity * orients[i].getY());
		positions[i][2] += dt * (forces[i][2] + activity * orients[i].getZ());
		positions[i][0] += noiseTemp(rng);
		positions[i][1] += noiseTemp(rng); 
		positions[i][2] += noiseTemp(rng); 
		// Rotational diffusion
		orients[i].randomRotation(stddevOrient, rng);
	}

	enforcePBC();
}

/* \brief Compute the forces between the particles.
 *
 * Implement harmonic spheres.
 */
void State3d::calcInternalForces() {
    for (long i = 0 ; i < n_parts ; ++i) {
		forces[i][0] = 0;
		forces[i][1] = 0;
		forces[i][2] = 0;
    }

	// Recompute the boxes
	boxes.update(&positions);
	const long n_boxes = boxes.getNBoxes();
	const std::vector< std::vector<long> > * nbrs_pos = boxes.getNbrsPos();
	const std::vector< std::vector<long> > * parts_of_box = \
		boxes.getPartsOfBox();

	for (long b1 = 0 ; b1 < n_boxes ; ++b1) {
		for (long b2 : (*nbrs_pos)[b1]) {
			for (long i : (*parts_of_box)[b1]) {
				for (long j : (*parts_of_box)[b2]) {
					double dx = positions[i][0] - positions[j][0];
					double dy = positions[i][1] - positions[j][1];
					double dz = positions[i][2] - positions[j][2];
					// We want the periodized interval to be centered in 0
					pbcSym(dx, len);
					pbcSym(dy, len);
					pbcSym(dz, len);
					double dr2 = dx * dx + dy * dy + dz * dz;

					if(dr2 * (1. - dr2) > 0.) {
						double u = pot_strength * (1.0 / std::sqrt(dr2) - 1.0);
						double fx = u * dx;
						double fy = u * dy;
						double fz = u * dz;

						forces[i][0] += fx;
						forces[j][0] -= fx;
						forces[i][1] += fy;
						forces[j][1] -= fy;
						forces[i][2] += fz;
						forces[j][2] -= fz;
					}
				}
			}
		}
	}
}

/* 
 * \brief Enforce periodic boundary conditions
 */
void State3d::enforcePBC() {
	for (long i = 0 ; i < n_parts ; ++i) {
		pbc(positions[i][0], len);
		pbc(positions[i][1], len);
		pbc(positions[i][2], len);
	}
}
