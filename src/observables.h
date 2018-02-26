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
 * \file observables.h
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Computation and export of the observables
 *
 * Header file for observables.cpp.
*/

#ifndef ACTIVEBROWNIAN_OBSERVABLES_H
#define ACTIVEBROWNIAN_OBSERVABLES_H

#include <vector>
#include "state.h"

class Observables {
	public:
		Observables(const double len_, const long n_parts_,
				    const double step_r_, const long n_div_angle_);
		//! Compute the observables for a given state
		void compute(const State *state);
		//! Export to hdf5
		void writeH5(const std::string fname, double rho, long n_parts,
	                 double pot_strength, double temperature, double rot_dif,
				     double activity, double dt, long n_iters) const;

	private:
		const double len; //!< Length of the box
		const long n_parts; //!< Number of particles 
		const double step_r; //!< Size of spatial division
		const long n_div_angle; //!< Number of divisions for angle
		const long n_div_r; //! Number of divisions in x
		const long n_div_tot; //! Total number of divisions

		std::vector<long long> correls; //!< Correlations
};

#endif // ACTIVEBROWNIAN_OBSERVABLES_H
