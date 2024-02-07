#include <cmath>
#include <chrono>
#include <algorithm>
#include <iostream>
#include "state.h"

/*Initialize the state of the system: particles randomly placed in a 2d box.*/
State::State(const double _Lx, const double _Ly, const long _n_parts,
	         const double _diam, const double _trans_dif,
			 const double _rot_dif, const double _activity, const double _dt,
			 const double _pot_strength, Infered &_infered) :
	lengths({_Lx, _Ly, 2 * M_PI}),
	n_parts(_n_parts),
	diam(_diam),
	activity(_activity),
	dt(_dt),
	pot_strength(_pot_strength),
	infered(_infered),
	boxes(_Lx, _Ly, _n_parts, _diam),
#ifdef USE_MKL
	stddevs({std::sqrt(2.0 * _trans_dif * dt),
			 std::sqrt(2.0 * _trans_dif * dt),
			 std::sqrt(2.0 * _rot_dif * dt)})
#else
	// We seed the RNG with the current time
	rng(std::chrono::system_clock::now().time_since_epoch().count()),
	// Gaussian noise from the translational diffusivity
	noiseTemp(0.0, std::sqrt(2.0 * _trans_dif * dt)),
	// Gaussian noise from the rotational diffusivity
	noiseAngle(0.0, std::sqrt(2.0 * _rot_dif * dt))
#endif
{

	for (int k = 0 ; k < 3 ; ++k) {
		positions[k].resize(n_parts);
		forces[k].assign(n_parts, 0);
	}

#ifdef USE_MKL
	for (int k = 0 ; k < 3 ; ++k) {
		aux[k].resize(n_parts);
	}

	vslNewStream(&stream, VSL_BRNG_SFMT19937,
			std::chrono::system_clock::now().time_since_epoch().count());
	for (int i = 0 ; i < 3 ; ++i) {
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
					 positions[i].data(), 0, lengths[i]);
	}
#else
    std::uniform_real_distribution<double> rndPosX(0, lengths[0]);
    std::uniform_real_distribution<double> rndPosY(0, lengths[1]);
    std::uniform_real_distribution<double> rndAngle(0, lengths[2]);

	for (long i = 0 ; i < n_parts ; ++i) {
		positions[0][i] = rndPosX(rng);
		positions[1][i] = rndPosY(rng);
		positions[2][i] = rndAngle(rng);
	}
#endif
}

/* Evolve the system according to coupled Langevin equations. */
void State::evolve() {
	for (int k = 0 ; k < 3 ; ++k)
		for (long i = 0 ; i < n_parts ; ++i)
			forces[k][i] = 0.;

	calcInternalForces();

#ifdef USE_MKL
	vdSinCos(n_parts, positions[2].data(), aux[1].data(), aux[0].data());
	// Activity and forces
	cblas_daxpy(n_parts, activity, aux[0].data(), 1, forces[0].data(), 1);
	cblas_daxpy(n_parts, activity, aux[1].data(), 1, forces[1].data(), 1);
	for (int k = 0 ; k < 3 ; ++k) {
		cblas_daxpy(n_parts, dt, forces[k].data(), 1, positions[k].data(), 1);
		// Diffusion (translational or rotational)
		if (stddevs[k] > 0.) {
			vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n_parts,
						  aux[k].data(), 0, stddevs[k]);
			vdAdd(n_parts, positions[k].data(), aux[k].data(),
				  positions[k].data());
		}
	}
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

/* Compute the forces between the particles. */
void State::calcInternalForces() {
    /*for (long i = 0 ; i < n_parts ; ++i)
		for (long j = 0 ; j < i ; ++j)
			calcSoftForceIJ(i, j);*/

	// Recompute the boxes
	boxes.update(positions[0], positions[1]);
	const long n_boxes = boxes.getNBoxes();
	const auto &nbrs_pos = boxes.getNbrsPos();
	const auto &parts_of_box = boxes.getPartsOfBox();

	// Repulsive forces (with boxes)
	for (long b1 = 0 ; b1 < n_boxes ; ++b1) {
		for (auto it_i = parts_of_box[b1].cbegin() ;
			 it_i != parts_of_box[b1].cend() ; ++it_i) {
			// Same box
			for (auto it_j = parts_of_box[b1].cbegin() ;
				 it_j != it_i ; ++it_j) {
				calcSoftForceIJ(*it_i, *it_j);
			}
			// Neighboring boxes
			for (long b2 : nbrs_pos[b1]) {
				for (auto it_j = parts_of_box[b2].cbegin() ;
					 it_j != parts_of_box[b2].cend() ; ++it_j) {
					calcSoftForceIJ(*it_i, *it_j);
				}
			}
		}
	}

	// Infered forces
    for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = 0 ; j < n_parts ; ++j) {
			if (i != j) {
				calcInferedForceIJ(i, j);
			}
		}
	}
}

//! Compute internal force between particles i and j (soft potential)
void State::calcSoftForceIJ(const long i, const long j) {
	double dx = positions[0][i] - positions[0][j];
	double dy = positions[1][i] - positions[1][j];

	// We want the periodized interval to be centered in 0
	pbcSym(dx, lengths[0]);
	pbcSym(dy, lengths[1]);
	double dr2 = (dx * dx + dy * dy) / (diam * diam);

	if (dr2 * (1. - dr2) > 0.) {
		double u = pot_strength * (1.0 / std::sqrt(dr2) - 1.0);
		double fx = u * dx;
		double fy = u * dy;

		forces[0][i] += fx;
		forces[0][j] -= fx;
		forces[1][i] += fy;
		forces[1][j] -= fy;
	}
}

void State::calcInferedForceIJ(const long i, const long j) {
	double dx = positions[0][i] - positions[0][j];
	double dy = positions[1][i] - positions[1][j];

	// We want the periodized interval to be centered on 0
	pbcSym(dx, lengths[0]);
	pbcSym(dy, lengths[1]);
/*#ifdef USE_MKL
	double dr, alpha;
	dr = dx * dx + dy * dy;
	vdSqrt(1, &dr, &dr);
	vdAtan2(1, &dy, &dx, &alpha); // Not that faster than std::atan2
#else*/
	double dr = std::sqrt(dx * dx + dy * dy);
	double alpha = std::atan2(dy, dx);
	double c = dx / dr;
	double s = dy / dr;
	double t1 = positions[2][i] - alpha;
	double t2 = positions[2][j] - alpha;

	double fr, ft, fo;
	infered.computeForces(dr, t1, t2, fr, ft, fo);

	forces[0][i] += fr * c - ft * s;
	forces[1][i] += fr * s + ft * c;
	forces[2][i] += fo;
}

/*
//! Compute internal force between particles i and j (WCA potential)
void State::calcWCAForceIJ(const long i, const long j) {
	double dx = positions[0][i] - positions[0][j];
	double dy = positions[1][i] - positions[1][j];
	// We want the periodized interval to be centered in 0
	pbcSym(dx, lengths[0]);
	pbcSym(dy, lengths[1]);
	double dr2 = (dx * dx + dy * dy) / (TWOONESIXTH * diam);

	if(dr2 * (TWOONESIXTH - dr2) > 0.) {
		double u = pot_strength * (48. * pow(dr2, -7.) - 24.*pow(dr2, -4.)); 
		double fx = u * dx;
		double fy = u * dy;

		forces[0][i] += fx;
		forces[0][j] -= fx;
		forces[1][i] += fy;
		forces[1][j] -= fy;
	}
}
*/

/* Enforce periodic boundary conditions */
void State::enforcePBC() {
	for (int k = 0 ; k < 3 ; ++k) {
#ifdef USE_MKL
		pbcMKL(positions[k], lengths[k], aux[k], n_parts);
#else
		for (long i = 0 ; i < n_parts ; ++i) {
			pbc(positions[k][i], lengths[k]);
		}
#endif
	}
}

/* Store positions and angles into 'out' */
void State::store(std::vector<double> &out) const {
	out.reserve(3 * n_parts);
	for (int i = 0 ; i < 3 ; ++i) {
		out.insert(out.end(), positions[i].begin(), positions[i].end());
	}
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
