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
#ifdef USE_MKL
	//aux.resize(n_modes_tot);
	aux.resize(n_funs);
#endif
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

	/*void cblas_dgemv (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans,
		const MKL_INT m, const MKL_INT n, const double alpha,
		const double *a, const MKL_INT lda, const double *x, const MKL_INT incx,
		const double beta, double *y, const MKL_INT incy); */

	/*
	// First solution, without transposition
	cblas_dgemv(CblasRowMajor, CblasNoTrans,
			    n_modes_tot, n_funs, 1.,
				coeffs[0].data(), n_funs, radialCompos.data(), 1,
				0., aux.data(), 1);
	f_r = cblas_ddot(n_modes_tot, aux.data(), 1, aux_cos.data(), 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans,
			    n_modes_tot, n_funs, 1.,
				coeffs[1].data(), n_funs, radialCompos.data(), 1,
				0., aux.data(), 1);
	f_t = cblas_ddot(n_modes_tot, aux.data(), 1, aux_sin.data(), 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans,
			    n_modes_tot, n_funs, 1.,
				coeffs[2].data(), n_funs, radialCompos.data(), 1,
				0., aux.data(), 1);
	f_o = cblas_ddot(n_modes_tot, aux.data(), 1, aux_sin.data(), 1);*/

	// Second solution, with transposition
	cblas_dgemv(CblasRowMajor, CblasTrans,
			    n_modes_tot, n_funs, 1.,
				coeffs[0].data(), n_funs, aux_cos.data(), 1,
				0., aux.data(), 1);
	f_r = cblas_ddot(n_funs, aux.data(), 1, radialCompos.data(), 1);
	cblas_dgemv(CblasRowMajor, CblasTrans,
			    n_modes_tot, n_funs, 1.,
				coeffs[1].data(), n_funs, aux_sin.data(), 1,
				0., aux.data(), 1);
	f_t = cblas_ddot(n_funs, aux.data(), 1, radialCompos.data(), 1);
	cblas_dgemv(CblasRowMajor, CblasTrans,
			    n_modes_tot, n_funs, 1.,
				coeffs[2].data(), n_funs, aux_sin.data(), 1,
				0., aux.data(), 1);
	f_o = cblas_ddot(n_funs, aux.data(), 1, radialCompos.data(), 1);

#else
	for (long i = 0 ; i < n_modes_tot ; ++i) {
		aux_cos[i] = std::cos(angles[i]);
		aux_sin[i] = std::sin(angles[i]);
	}

	// Forces
	f_r = f_t = f_o = 0.;
	for (long j = 0 ; j < n_modes_tot ; ++j) {
		for (long k = 0 ; k < n_funs ; ++k) {
			f_r += coeffs[0][j*n_funs+k] * radialCompos[k] * aux_cos[j];
			f_t += coeffs[1][j*n_funs+k] * radialCompos[k] * aux_sin[j];
			f_o += coeffs[2][j*n_funs+k] * radialCompos[k] * aux_sin[j];
		}
	}
#endif
}

// Store (r/r0)^k exp(-r/r0) / k! for 0 <= k < n_funs
void Infered::computeRadialCompos(double r) {
	r /= r0;
	radialCompos[0] = std::exp(-r);
	for (int i = 1 ; i < n_funs ; ++i) {
		radialCompos[i] = radialCompos[i-1] * r / i;
	}
}
