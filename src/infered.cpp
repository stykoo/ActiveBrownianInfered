#include <iostream>
#include <array>
#include <vector>
#include <cassert>
#include "infered.h"

#ifdef USE_MKL
	#include "mkl.h"
#endif

Infered::Infered(const std::array<std::vector<int>, 2> &_ks,
		         const std::array<std::vector<double>, 3> &_coeffs,
				 double _r0) :
		ks(_ks), coeffs(_coeffs), n_modes_tot(_ks[0].size()),
		n_funs(_coeffs[0].size() / _ks[0].size()), r0(_r0) {
	radialCompos.resize(n_funs);
	angles.resize(n_modes_tot);
	aux_cos.resize(n_modes_tot);
	aux_sin.resize(n_modes_tot);
}

void Infered::computeForces(const double r, const double t1, const double t2,
		                    double &f_r, double &f_t, double &f_o) {
	// Radial components
	computeRadialCompos(r);

	// Angular components
	for (long i = 0 ; i < n_modes_tot ; ++i) {
		angles[i] = ks[0][i] * t1 + ks[1][i] * t2; // angles
	}

	// Sin and cos of angles
#ifdef USE_MKL
	vdSinCos(n_modes_tot, angles.data(), aux_sin.data(), aux_cos.data());
#else
	for (long i = 0 ; i < n_modes_tot ; ++i) {
		aux_cos[i] = std::cos(angles[i]);
		aux_sin[i] = std::sin(angles[i]);
	}
#endif

	// Forces
	f_r = f_t = f_o = 0.;
	for (long j = 0 ; j < n_modes_tot ; ++j) {
		for (long k = 0 ; k < n_funs ; ++k) {
			f_r += coeffs[0][j*n_funs+k] * radialCompos[k] * aux_cos[j];
			f_t += coeffs[1][j*n_funs+k] * radialCompos[k] * aux_sin[j];
			f_o += coeffs[2][j*n_funs+k] * radialCompos[k] * aux_sin[j];
		}
	}
}

// Store (r/r0)^k exp(-r/r0) / k! for 0 <= k < n_funs
void Infered::computeRadialCompos(const double r) {
	radialCompos[0] = std::exp(-r / r0);
	for (int i = 1 ; i < n_funs ; ++i) {
		radialCompos[i] = radialCompos[i-1] * (r / r0) / i;
	}
}
