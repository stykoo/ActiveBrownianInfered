#ifndef ACTIVEBROWNIAN_INFERED_H_
#define ACTIVEBROWNIAN_INFERED_H_

#include <string>
#include <vector>
#include <cmath>

class Infered {
	public:
		Infered(long n_funs_, long n_modes_, std::string fname_r,
				std::string fname_t, std::string fname_o);
		void computeForces(const double r, const double t1, const double t2,
		                   double &f_r, double &f_t, double &f_o);

	private:
		void loadCoeffs(std::string fname,
				        std::vector<std::vector<double>> &coeffs);
		void computeRadialCompos(const double r);
		void computeAngularCompos(const double t1, const double t2);

		const long n_funs; //!< Number of base functions
		const long n_modes; //!< Number of angular modes
		const long n_modes_tot; //!< Total number of angular modes

		std::vector< std::vector<double> > coeffs_r; //!< Coefficients
		std::vector< std::vector<double> > coeffs_t; //!< Coefficients
		std::vector< std::vector<double> > coeffs_o; //!< Coefficients
		//std::vector<double> coeffs_r; //!< Coefficients
		//std::vector<double> coeffs_t; //!< Coefficients
		//std::vector<double> coeffs_o; //!< Coefficients

		std::vector<long> ks1; //!< Coefficients in front of t1 in cos
		std::vector<long> ks2; //!< Coefficients in front of t2 in cos
		std::vector<double> radialCompos; //!< Radial components
		std::vector<double> angles;
		std::vector<double> aux_cos;
		std::vector<double> aux_sin;
};

//void base_funs(const double r, const int n_funs, std::vector<double> &res);

#endif
