#include <fstream>
#include <iterator>
#include <sstream>
#include <cassert>
#include "infered.h"

#ifdef USE_MKL
	#include "mkl.h"
#endif

Infered::Infered(long n_funs_, long n_modes_, std::string fname_r,
	std::string fname_t, std::string fname_o) :
		n_funs(n_funs_), n_modes(n_modes_),
		n_modes_tot(n_modes_ * (n_modes_ + 1)) { // Load coefficients
	loadCoeffs(fname_r, coeffs_r);
	loadCoeffs(fname_t, coeffs_t);
	loadCoeffs(fname_o, coeffs_o);

	radialCompos.resize(n_funs);
	angles.resize(n_modes_tot);
	aux_cos.resize(n_modes_tot);
	aux_sin.resize(n_modes_tot);

	// Compute coefficients in front of theta1 and theta2 in cos
	for (long n = 0 ; n < n_modes ; ++n) {
		for (long k = -n+1 ; k <=n ; ++k) {
			ks1.push_back(n-std::abs(k));
			ks2.push_back(k);
		}
	}
	assert(ks1.size() == (size_t) n_modes_tot);
	assert(ks2.size() == (size_t) n_modes_tot);
}

void Infered::loadCoeffs(std::string fname,
		                 std::vector<std::vector<double>> &coeffs) {
	// TODO: move to std::vector<double> (1D)
	
	std::string line;
	std::ifstream infile(fname);
	assert(infile);

	// Read matrix of coefficients
	while (std::getline(infile, line)) {
        std::istringstream ss(line);
        coeffs.emplace_back(
				std::istream_iterator<double>(ss),
				std::istream_iterator<double>()
			);
    }

	// Check if dimensions are correct
	assert(coeffs.size() == (size_t) n_funs);
	for (long i = 0 ; i < n_funs ; ++i) {
		assert(coeffs[i].size() == (size_t) n_modes_tot);
	}
}

void Infered::computeForces(const double r, const double t1, const double t2,
		                    double &f_r, double &f_t, double &f_o) {
	// Radial components
	computeRadialCompos(r);

	// Angular components
	for (long i = 0 ; i < n_modes_tot ; ++i) {
		angles[i] = ks1[i] * t1 + ks2[i] * t2; // angles
	}
#ifdef USE_MKL
	vdSinCos(n_modes_tot, angles.data(), aux_sin.data(), aux_cos.data());
#else
	for (long i = 0 ; i < n_modes_tot ; ++i) {
		aux_cos[i] = std::cos(angles[i]); // cos(angles)
		aux_sin[i] = std::sin(angles[i]); // sin(angles)
	}
#endif

/*#ifdef USE_MKL
	// TODO: use cblas_dgemv and cblas_ddot
	cblas_dgemv(CblasColMajor, CblasNoTrans, n_modes_tot, n_funs,
			    1., coeffs_r.data(), n_modes_tot, radialCompos.data(), 1,
				0., angles.data(), 1); // Use angles as temporary array
	f_r = cblas_ddot(n_modes_tot, aux_cos.data(), 1, angles.data(), 1);
#else */
	f_r = f_t = f_o = 0.;
	for (long i = 0 ; i < n_funs ; ++i) {
		for (long j = 0 ; j < n_modes_tot ; ++j) {
			f_r += coeffs_r[i][j] * radialCompos[i] * aux_cos[j];
			f_t += coeffs_t[i][j] * radialCompos[i] * aux_sin[j];
			f_o += coeffs_o[i][j] * radialCompos[i] * aux_sin[j];
			//f_r += coeffs_r[i*nmodes_tot+j] * radialCompos[i] * aux_cos[j];
			//f_t += coeffs_t[i*nmodes_tot+j] * radialCompos[i] * aux_sin[j];
			//f_o += coeffs_o[i*nmodes_tot+j] * radialCompos[i] * aux_sin[j];
		}
	}
// #endif
}

// Store r^k exp(-r) / k! for 0 <= k < n_funs
void Infered::computeRadialCompos(const double r) {
	radialCompos[0] = std::exp(-r);
	for (int i = 1 ; i < n_funs ; ++i) {
		radialCompos[i] = radialCompos[i-1] * r / i;
	}
}
