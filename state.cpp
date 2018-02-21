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

/*
 * \brief Constructor of boxes.
 *
 * Construct the boxes with which we will be able to classify the particles.
 *
 * \param len Length of the system
 * \param n_parts Number of particles
 */
Boxes::Boxes(const double len, const long n_parts) :
		n_boxes_x(std::floor(len)), n_boxes(n_boxes_x * n_boxes_x),
		len_box(len / n_boxes_x), n_parts(n_parts) {
	nbrs.resize(n_boxes);
	computeNbrs();
	box_of_part.resize(n_parts);
	parts_of_box.resize(n_boxes);
}

/*
 * \brief Classify the particles of a state in the boxes.
 * 
 * Warning: the state should be consistent with the Boxes object!
 *
 * \param state State of the system
 */
void Boxes::update(const State *state) {
    for (long i=0 ; i < n_boxes ; ++i) {
		parts_of_box[i].clear();
	}

    for (long i=0 ; i < n_parts ; ++i) {
		long bx = (long) std::floor(state->getPosX(i) / len_box);
		long by = (long) std::floor(state->getPosY(i) / len_box);
		long box = by * n_boxes_x + bx;
		box_of_part[i] = box;
		parts_of_box[box].push_front(i);
	}
}

/*
 * \brief Compute the numbers of the neighboring boxes of each box.
 */
void Boxes::computeNbrs() {
	for (long i = 0 ; i < n_boxes ; ++i) {
		long xs[3];
		long ys[3];
		xs[1] = i % n_boxes_x;
		ys[1] = i / n_boxes_x;
		xs[0] = (xs[1] + n_boxes_x - 1) % n_boxes_x;
		ys[0] = (ys[1] + n_boxes_x - 1) % n_boxes_x;
		xs[2] = (xs[1] + 1) % n_boxes_x;
		ys[2] = (ys[1] + 1) % n_boxes_x;

		for (int k = 0 ; k < 9 ; ++k) {
			nbrs[i][k] = ys[k / 3] * n_boxes_x + xs[k % 3];
		}
	}
}

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
	calcInternalForces();

	for (long i = 0 ; i < n_parts ; ++i) {
		// Internal forces +  Activity + Gaussian noise
		positions[i][0] += dt * (forces[i][0]
		                         + activity * std::cos(angles[i]));
		positions[i][1] += dt * (forces[i][1]
		                         + activity * std::sin(angles[i]));
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
	boxes.update(this);

    for (long i = 0 ; i < n_parts ; ++i) {
		long k = boxes.getBoxOfPart(i);
		/*std::cout << "(" << positions[i][0] << ", " << positions[i][1]
			      << ") -> " << k << "\n";*/

		for(auto it_b = boxes.getNbrsBegin(k) ; it_b != boxes.getNbrsEnd(k) ; 
		    ++it_b) {
            for (auto it_j = boxes.getPartsOfBoxBegin(*it_b) ;
                 it_j != boxes.getPartsOfBoxEnd(*it_b) && (*it_j) > i ;
                 ++it_j) {
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
