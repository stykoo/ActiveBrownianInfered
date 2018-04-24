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
#include <algorithm>
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
	boxes(_len, _n_parts),
#ifdef USE_MKL
	stddev_temp(std::sqrt(2.0 * _temperature * dt)),
	stddev_rot(std::sqrt(2.0 * _rot_dif * dt))
#else
	// We seed the RNG with the current time
	rng(std::chrono::system_clock::now().time_since_epoch().count()),
	// Gaussian noise from the temperature
	noiseTemp(0.0, std::sqrt(2.0 * _temperature * dt)),
	// Gaussian noise from the rotational diffusivity
	noiseAngle(0.0, std::sqrt(2.0 * _rot_dif * dt))
#endif
{
	positions.resize(n_parts);
	angles.resize(n_parts);
	forces.resize(2 * n_parts);

#ifdef USE_MKL
	vslNewStream(&stream, VSL_BRNG_SFMT19937,
			std::chrono::system_clock::now().time_since_epoch().count());
	ran_num_x.resize(n_parts);
	ran_num_y.resize(n_parts);
	ran_num_angle.resize(n_parts);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts, ran_num_x.data(),
			     0, len);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts, ran_num_y.data(),
			     0, len);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
			     ran_num_angle.data(), 0, 2.0 * M_PI);
#else
    std::uniform_real_distribution<double> rndPos(0, len);
    std::uniform_real_distribution<double> rndAngle(0, 2.0 * M_PI);
#endif

	for (long i = 0 ; i < n_parts ; ++i) {
#ifdef USE_MKL
		positions[i][0] = ran_num_x[i];
		positions[i][1] = ran_num_y[i];
		angles[i] = ran_num_angle[i];
#else
		positions[i][0] = rndPos(rng);
		positions[i][1] = rndPos(rng);
		angles[i] = rndAngle(rng);
#endif
		forces[i] = 0;
		forces[n_parts + i] = 0;
	}
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation.
 */
void State::evolve() {
	double c, s;
	calcInternalForces();
#ifdef USE_MKL
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n_parts,
			      ran_num_x.data(), 0, stddev_temp);
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n_parts,
			      ran_num_y.data(), 0, stddev_temp);
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n_parts,
			      ran_num_angle.data(), 0, stddev_rot);
#endif

	for (long i = 0 ; i < n_parts ; ++i) {
		// Computation of sin and cos
#ifdef __GNUC__
		sincos(angles[i], &s, &c);
#else
		s = std::sin(angles[i]);
		c = std::cos(angles[i]);
#endif
		// Internal forces +  Activity + Gaussian noise
		positions[i][0] += dt * (forces[i] + activity * c);
		positions[i][1] += dt * (forces[n_parts + i] + activity * s);
		// Diffusion and rotational diffusion
#ifdef USE_MKL
		positions[i][0] += ran_num_x[i];
		positions[i][1] += ran_num_y[i];
		angles[i] += ran_num_angle[i];
#else
		positions[i][0] += noiseTemp(rng);
		positions[i][1] += noiseTemp(rng); 
		angles[i] += noiseAngle(rng);
#endif
	}

	enforcePBC();
}

/* \brief Compute the forces between the particles.
 *
 * Implement harmonic spheres.
 */
void State::calcInternalForces() {
    for (long i = 0 ; i < 2 * n_parts ; ++i) {
		forces[i] = 0;
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
					// We want the periodized interval to be centered in 0
					pbcSym(dx, len);
					pbcSym(dy, len);
					double dr2 = dx * dx + dy * dy;

					if(dr2 * (1. - dr2) > 0.) {
						double u = pot_strength * (1.0 / std::sqrt(dr2) - 1.0);
						double fx = u * dx;
						double fy = u * dy;

						forces[i] += fx;
						forces[j] -= fx;
						forces[n_parts + i] += fy;
						forces[n_parts + j] -= fy;
					}
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
