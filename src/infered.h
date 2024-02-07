#ifndef ACTIVEBROWNIAN_INFERED_H_
#define ACTIVEBROWNIAN_INFERED_H_

#include <string>
#include <vector>
#include <cmath>

class Infered {
	public:
		Infered(const std::array<std::vector<int>, 2> &_ks,
		  	    const std::array<std::vector<double>, 3> &_coeffs,
				double _r0);
		void computeForces(const double r, const double t1, const double t2,
		                   double &f_r, double &f_t, double &f_o);

	private:
		void computeRadialCompos(double r);

		const std::array<std::vector<int>, 2> &ks; //!< Factors for cos/sin
		const std::array<std::vector<double>, 3> &coeffs;  //!< Coefficients
		const long n_modes_tot; //!< Total number of angular modes
		const long n_funs; //!< Number of base functions
		const double r0; //!< Radial increment

		std::vector<double> radialCompos; //!< Radial components
		std::vector<double> angles;
		std::vector<double> aux_cos;
		std::vector<double> aux_sin;
#ifdef USE_MKL
		std::vector<double> aux;
#endif
};

#endif
