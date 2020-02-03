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
#include <iostream>
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
 * \param _fac_boxes Factor for the boxes
 */
State::State(const double _len, const long _n_parts,
	         const double _pot_strength, const double _temperature,
			 const double _rot_dif, const double _activity, const double _dt,
			 const int _fac_boxes, const bool _wca) :
	len(_len), n_parts(_n_parts), pot_strength(_pot_strength),
	activity(_activity), dt(_dt), wca(_wca),
	boxes(_len, _n_parts, (_wca ? TWOONESIXTH : 1.0), _fac_boxes),
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
	positions[0].resize(n_parts);
	positions[1].resize(n_parts);
	angles.resize(n_parts);
	forces[0].assign(n_parts, 0);
	forces[1].assign(n_parts, 0);
	f_along.assign(n_parts, 0);

#ifdef USE_MKL
	aux_x.resize(n_parts);
	aux_y.resize(n_parts);
	aux_angle.resize(n_parts);

	vslNewStream(&stream, VSL_BRNG_SFMT19937,
			std::chrono::system_clock::now().time_since_epoch().count());
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
			     positions[0].data(), 0, len);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
			     positions[1].data(), 0, len);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
			     angles.data(), 0, 2.0 * M_PI);
#else
    std::uniform_real_distribution<double> rndPos(0, len);
    std::uniform_real_distribution<double> rndAngle(0, 2.0 * M_PI);

	for (long i = 0 ; i < n_parts ; ++i) {
		positions[0][i] = rndPos(rng);
		positions[1][i] = rndPos(rng);
		angles[i] = rndAngle(rng);
		forces[0][i] = 0;
		forces[1][i] = 0;
	}
#endif
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation.
 */
void State::evolve() {
	calcInternalForces();

#ifdef USE_MKL
	vdSinCos(n_parts, angles.data(), aux_y.data(), aux_x.data());
	// It doesn't look optimal to do it each time
	// but the alternative would be to store a copy of the angles.
	for (long i = 0 ; i < n_parts ; ++i) {
		f_along[i] = forces[0][i] * aux_x[i] + forces[1][i] * aux_y[i];
	}
	// Activity and forces
	cblas_daxpy(n_parts, activity, aux_x.data(), 1, forces[0].data(), 1);
	cblas_daxpy(n_parts, activity, aux_y.data(), 1, forces[1].data(), 1);
	cblas_daxpy(n_parts, dt, forces[0].data(), 1, positions[0].data(), 1);
	cblas_daxpy(n_parts, dt, forces[1].data(), 1, positions[1].data(), 1);
	// Diffusion and rotational diffusion
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n_parts,
			      aux_x.data(), 0, stddev_temp);
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n_parts,
			      aux_y.data(), 0, stddev_temp);
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n_parts,
			      aux_angle.data(), 0, stddev_rot);
	vdAdd(n_parts, positions[0].data(), aux_x.data(), positions[0].data());
	vdAdd(n_parts, positions[1].data(), aux_y.data(), positions[1].data());
	vdAdd(n_parts, angles.data(), aux_angle.data(), angles.data());
#else
	double c, s;
	for (long i = 0 ; i < n_parts ; ++i) {
		// Computation of sin and cos
	#ifdef __GNUC__
		sincos(angles[i], &s, &c);
	#else
		s = std::sin(angles[i]);
		c = std::cos(angles[i]);
	#endif
		f_along[i] = forces[0][i] * c + forces[1][i] * s;
		// Internal forces +  Activity + Gaussian noise
		positions[0][i] += dt * (forces[0][i] + activity * c);
		positions[1][i] += dt * (forces[1][i] + activity * s);
		// Diffusion and rotational diffusion
		positions[0][i] += noiseTemp(rng);
		positions[1][i] += noiseTemp(rng); 
		angles[i] += noiseAngle(rng);
	}
#endif

	enforcePBC();
}

double State::avgFAlong() const {
	double f = 0.0;
	for (long i = 0 ; i < n_parts ; ++i) {
		f += f_along[i];
	}
	return f / n_parts;
}

void State::dump() const {
	for (long i = 0 ; i < n_parts ; ++i) {
		std::cout << positions[0][i] << " "
			<< positions[1][i] << " "
			<< angles[i] * 180 / M_PI << "\n";
	}
}

/* \brief Compute the forces between the particles.
 *
 * Implement harmonic spheres.
 */
void State::calcInternalForces() {
    for (long i = 0 ; i < n_parts ; ++i) {
		forces[0][i] = 0;
		forces[1][i] = 0;
    }

	// Recompute the boxes
	boxes.update(positions);
	const long n_boxes = boxes.getNBoxes();
	const std::vector< std::vector<long> > &nbrs_pos = boxes.getNbrsPos();
	const std::vector< std::vector<long> > &parts_of_box = \
		boxes.getPartsOfBox();

	/*for (long b1 = 0 ; b1 < n_boxes ; ++b1) {
		std::cout << "-> " << b1 << ": ";
		for (auto it_i = parts_of_box[b1].cbegin() ;
			 it_i != parts_of_box[b1].cend() ; ++it_i) {
			std::cout << *it_i << " ";
		}
		std::cout << std::endl;
	}*/

	if (wca) {
		for (long b1 = 0 ; b1 < n_boxes ; ++b1) {
			//std::cout << "-> " << b1 << "\n";
			for (auto it_i = parts_of_box[b1].cbegin() ;
				 it_i != parts_of_box[b1].cend() ; ++it_i) {
				// Same box
				for (auto it_j = parts_of_box[b1].cbegin() ;
					 it_j != it_i ; ++it_j) {
					calcInternalForceIJ_WCA(*it_i, *it_j);
				}
				// Neighboring boxes
				for (long b2 : nbrs_pos[b1]) {
					//std::cout << "[" << b1 << ", " << b2 << "]\n";
					for (auto it_j = parts_of_box[b2].cbegin() ;
						 it_j != parts_of_box[b2].cend() ; ++it_j) {
						calcInternalForceIJ_WCA(*it_i, *it_j);
					}
				}
			}
		}
	} else {
		for (long b1 = 0 ; b1 < n_boxes ; ++b1) {
			for (auto it_i = parts_of_box[b1].cbegin() ;
				 it_i != parts_of_box[b1].cend() ; ++it_i) {
				// Same box
				for (auto it_j = parts_of_box[b1].cbegin() ;
					 it_j != it_i ; ++it_j) {
					calcInternalForceIJ_soft(*it_i, *it_j);
				}
				// Neighboring boxes
				for (long b2 : nbrs_pos[b1]) {
					for (auto it_j = parts_of_box[b2].cbegin() ;
						 it_j != parts_of_box[b2].cend() ; ++it_j) {
						calcInternalForceIJ_soft(*it_i, *it_j);
					}
				}
			}
		}
	}
}

//! Compute internal force between particles i and j (soft potential)
void State::calcInternalForceIJ_soft(const long i, const long j) {
	// std::cout << i << " " << j << "\n";

	double dx = positions[0][i] - positions[0][j];
	double dy = positions[1][i] - positions[1][j];
	// We want the periodized interval to be centered in 0
	pbcSym(dx, len);
	pbcSym(dy, len);
	double dr2 = dx * dx + dy * dy;

	if(dr2 * (1. - dr2) > 0.) {
		double u = pot_strength * (1.0 / std::sqrt(dr2) - 1.0);
		double fx = u * dx;
		double fy = u * dy;

		forces[0][i] += fx;
		forces[0][j] -= fx;
		forces[1][i] += fy;
		forces[1][j] -= fy;
	}
}

//! Compute internal force between particles i and j (WCA potential)
void State::calcInternalForceIJ_WCA(const long i, const long j) {
	//std::cout << i << " " << j << "\n";

	double dx = positions[0][i] - positions[0][j];
	double dy = positions[1][i] - positions[1][j];
	// We want the periodized interval to be centered in 0
	pbcSym(dx, len);
	pbcSym(dy, len);
	double dr2 = dx * dx + dy * dy;

	if(dr2 * (TWOONESIXTH - dr2) > 0.) {
		double u = pot_strength;
		u *= (48. * pow(dr2, -7.) - 24.*pow(dr2, -4.)); 
		double fx = u * dx;
		double fy = u * dy;

		forces[0][i] += fx;
		forces[0][j] -= fx;
		forces[1][i] += fy;
		forces[1][j] -= fy;
	}
}

/* 
 * \brief Enforce periodic boundary conditions
 */
void State::enforcePBC() {
#ifdef USE_MKL
	pbcMKL(positions[0], len, aux_x, n_parts);
	pbcMKL(positions[1], len, aux_y, n_parts);
	pbcMKL(angles, 2.0 * M_PI, aux_angle, n_parts);
#else
	for (long i = 0 ; i < n_parts ; ++i) {
		pbc(positions[0][i], len);
		pbc(positions[1][i], len);
		pbc(angles[i], 2.0 * M_PI);
	}
#endif
}

#ifdef USE_MKL
void pbcMKL(std::vector<double> &v, const double L, std::vector<double> &aux,
	        const long N) {
	cblas_daxpby(N, 1.0 / L, v.data(), 1, 0.0, aux.data(), 1);
	vdFloor(N, aux.data(), aux.data());
	cblas_daxpy(N, -L, aux.data(), 1, v.data(), 1);
}	

void pbcSymMKL(std::vector<double> &v, const double L,
		       std::vector<double> &aux, const long N) {
	cblas_daxpby(N, 1.0 / L, v.data(), 1, 0.0, aux.data(), 1);
	vdRound(N, aux.data(), aux.data());
	cblas_daxpy(N, -L, aux.data(), 1, v.data(), 1);
}	
#endif
