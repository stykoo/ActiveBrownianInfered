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
 * \param _temperature Temperature
 * \param _rot_dif Rotational diffusivity
 * \param _activity Activity
 * \param _dt Timestep
 * \param _fac_boxes Factor for the boxes
 */
State::State(const double _len, const long _n_parts,
	         const double _temperature,
			 const double _rot_dif, const double _activity, const double _dt,
			 Infered &_infered) :
	len(_len), n_parts(_n_parts),
	activity(_activity), dt(_dt), infered(_infered),
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
	for (int k = 0 ; k < 3 ; ++k) {
		positions[k].resize(n_parts);
		forces[k].assign(n_parts, 0);
	}

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
			     positions[2].data(), 0, 2.0 * M_PI);
#else
    std::uniform_real_distribution<double> rndPos(0, len);
    std::uniform_real_distribution<double> rndAngle(0, 2.0 * M_PI);

	for (long i = 0 ; i < n_parts ; ++i) {
		positions[0][i] = rndPos(rng);
		positions[1][i] = rndPos(rng);
		positions[2][i] = rndAngle(rng);
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
	vdSinCos(n_parts, positions[2].data(), aux_y.data(), aux_x.data());
	// Activity and forces
	cblas_daxpy(n_parts, activity, aux_x.data(), 1, forces[0].data(), 1);
	cblas_daxpy(n_parts, activity, aux_y.data(), 1, forces[1].data(), 1);
	cblas_daxpy(n_parts, dt, forces[0].data(), 1, positions[0].data(), 1);
	cblas_daxpy(n_parts, dt, forces[1].data(), 1, positions[1].data(), 1);
	cblas_daxpy(n_parts, dt, forces[2].data(), 1, positions[2].data(), 1);
	// Diffusion and rotational diffusion
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n_parts,
			      aux_x.data(), 0, stddev_temp);
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n_parts,
			      aux_y.data(), 0, stddev_temp);
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n_parts,
			      aux_angle.data(), 0, stddev_rot);
	vdAdd(n_parts, positions[0].data(), aux_x.data(), positions[0].data());
	vdAdd(n_parts, positions[1].data(), aux_y.data(), positions[1].data());
	vdAdd(n_parts, positions[2].data(), aux_angle.data(), positions[2].data());
#else
	double c, s;
	for (long i = 0 ; i < n_parts ; ++i) {
		// Computation of sin and cos
	#ifdef __GNUC__
		sincos(positions[2][i], &s, &c);
	#else
		s = std::sin(positions[2][i]);
		c = std::cos(positions[2][i]);
	#endif
		// Internal forces +  Activity + Gaussian noise
		positions[0][i] += dt * (forces[0][i] + activity * c);
		positions[1][i] += dt * (forces[1][i] + activity * s);
		positions[2][i] += dt * forces[2][i];
		// Diffusion and rotational diffusion
		positions[0][i] += noiseTemp(rng);
		positions[1][i] += noiseTemp(rng); 
		positions[2][i] += noiseAngle(rng);
	}
#endif

	enforcePBC();
}

/* \brief Compute the forces between the particles.
 *
 * Implement harmonic spheres.
 */
void State::calcInternalForces() {
    for (long i = 0 ; i < n_parts ; ++i) {
		forces[0][i] = 0;
		forces[1][i] = 0;
		forces[2][i] = 0;
    }

    for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = 0 ; j < n_parts ; ++j) {
			if (i != j) {
				calcInferedForceIJ(i, j);
			}
		}
	}
}

void State::calcInferedForceIJ(const long i, const long j) {
	double dx = positions[0][i] - positions[0][j];
	double dy = positions[1][i] - positions[1][j];

	// We want the periodized interval to be centered on 0
	pbcSym(dx, len);
	pbcSym(dy, len);
	double dr = std::sqrt(dx * dx + dy * dy);
	double c = dx / dr;
	double s = dy / dr;

	double fr, ft, fo;
	infered.computeForces(dr, positions[2][i], positions[2][j], fr, ft, fo);

	forces[0][i] += fr * c - ft * s;
	forces[1][i] += fr * s + ft * c;
	forces[2][i] += fo;
}

/* 
 * \brief Enforce periodic boundary conditions
 */
void State::enforcePBC() {
#ifdef USE_MKL
	pbcMKL(positions[0], len, aux_x, n_parts);
	pbcMKL(positions[1], len, aux_y, n_parts);
	pbcMKL(positions[2], 2.0 * M_PI, aux_angle, n_parts);
#else
	for (long i = 0 ; i < n_parts ; ++i) {
		pbc(positions[0][i], len);
		pbc(positions[1][i], len);
		pbc(positions[2][i], 2.0 * M_PI);
	}
#endif
}

void State::dump() const {
	for (long i = 0 ; i < n_parts ; ++i) {
		std::cout << positions[0][i] << " "
			<< positions[1][i] << " "
			<< positions[2][i] * 180 / M_PI << "\n";
	}
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
