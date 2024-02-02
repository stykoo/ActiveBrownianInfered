#ifndef ACTIVEBROWNIAN_OBSERVABLES_H
#define ACTIVEBROWNIAN_OBSERVABLES_H

#include <vector>
#include "state.h"
#include "infered.h"

class Observables {
	public:
		Observables(const double len_, const long n_parts_,
		            const double step_r, const long n_div_angle_);
		//! Compute the observables for a given state
		void compute(const State *state);
		//! Compute the forces on the grid
		void computeForces(Infered &infered);
		//! Export to hdf5
		void writeH5(const std::string fname, double rho, long n_parts,
	                 double temperature, double rot_dif,
				     double activity, double dt, long n_iters, long n_iters_th,
					 long skip) const;

	private:
		const double len; //!< Length of the box
		const long n_parts; //!< Number of particles 
		double step_r; //!< Size of spatial division
		const long n_div_angle; //!< Number of divisions for angle
		double scal_r; //!< Scale for spatial divisions
		const double step_angle; //!< Step for angular divisions
		long n_div_r; //!< Number of divisions in x
		long n_div_tot; //!< Total number of divisions
		
		std::vector<long long> forces_r, forces_t, forces_o;

#ifdef USE_MKL
		const long n_pairs;
		//const std::vector<double> ones;
		std::vector<double> dxs, dys, phis, drs, thetas1;
#endif

		long n_calls; //!< Number of calls of 'compute'
		std::vector<long long> correls; //!< Correlations
		std::array<std::vector<double>, 3> forces; 
};

#endif // ACTIVEBROWNIAN_OBSERVABLES_H
