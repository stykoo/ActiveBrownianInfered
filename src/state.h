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
 * \file state.h
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief State of the system
 *
 * Header file for state.h.
 * It defines the class State.
 */

#ifndef ACTIVEBROWNIAN_STATE_H_
#define ACTIVEBROWNIAN_STATE_H_

#include <vector>
#include <array>
#include "infered.h"
//#include "boxes.h"

#ifdef USE_MKL
	#include "mkl.h"
	#include "mkl_vsl.h"
#else
	#include <random>
#endif

// 2^(1/6)
#define TWOONESIXTH 1.12246204830937298143 


/*!
 * \brief Class for the state of the system
 *
 * This class takes care of the initialization
 * and the evolution of the state of the system.
 */
class State {
	public:
		//! Constructor of State
		State(const double _len, const long _n_parts, const double _temperature,
			  const double _rot_dif, const double _activity, const double _dt,
			  Infered &_infered);
		~State() {
#ifdef USE_MKL
			vslDeleteStream(&stream);
#endif
		}
		void evolve(); //!< Do one time step

		//! Get the x coordinate of the positions 
		const std::vector<double> & getPosX() const {
			return positions[0];
		}
		//! Get the y coordinate of the positions 
		const std::vector<double> & getPosY() const {
			return positions[1];
		}
		//! Get angle of particle i
		const std::vector<double> & getAngles() const {
			return angles;
		}

		void dump() const; //!< Dump the positions and orientations


	private:
		void calcInternalForces(); //!< Compute internal forces
		void calcInferedForceIJ(const long i, const long j);
		void enforcePBC(); //!< Enforce periodic boundary conditions

		const double len; //!< Length of the box
		const long n_parts; //!< Number of particles
		const double activity; //!< Activity
		const double dt; //!< Timestep
		Infered &infered; //!< Structure for inference

#ifdef USE_MKL
		double stddev_temp, stddev_rot;
		VSLStreamStatePtr stream;
		std::vector<double> aux_x, aux_y, aux_angle;
#else
		std::mt19937 rng; //!< Random number generator
		//! Gaussian noise for temperature
		std::normal_distribution<double> noiseTemp;
		//! Gaussian noise for angle
		std::normal_distribution<double> noiseAngle;
#endif

		//! Positions of the particles
		std::array<std::vector<double>, 2> positions;
		std::vector<double> angles; //<! Angles
		std::array<std::vector<double>, 3> forces;  //!< Internal forces
};

/*! 
 * \brief Periodic boundary conditions on a segment
 * 
 * Update x to be between 0 and L.
 *
 * \param x Value
 * \param L Length of the box
 */
template<typename T>
void pbc(T &x, const T L){
	x -= L * std::floor(x / L);
}
/*template<>
void pbc<double>(double &x, const double L) {
	// Trick to avoid floor
	double a = x / L;
	long i = (long) a;
	x -= L * (i - (i > a));
}*/

#ifdef USE_MKL
void pbcMKL(std::vector<double> &v, const double L, std::vector<double> &aux,
	        const long N);
void pbcSymMKL(std::vector<double> &v, const double L,
		       std::vector<double> &aux, const long N);
#endif

/*! 
 * \brief Periodic boundary conditions on a segment (symmetric)
 * 
 * Update x to be between -L/2 and L/2.
 *
 * \param x Value
 * \param L Length of the box
 */
template<typename T>
inline void pbcSym(T &x, const T L) {
	x -= L * std::round(x / L);
}

// Trick to avoid round ASSUMING LITTLE ENDIAN
union i_cast {double d; int i[2];};
#define double2int(i, d, t)  \
    {volatile union i_cast u; u.d = (d) + 6755399441055744.0; \
    (i) = (t)u.i[0];}

template<>
inline void pbcSym<double>(double &x, const double L) {
	double d = x / L;
	int i;
	double2int(i, d, int);
	x -= L * i;
}

/*! 
 * \brief Periodic boundary conditions on a segment, with offset
 * 
 * Update x to be between 0 and L, and o to be the corresponding offset
 *
 * \param x Value
 * \param o Offset
 * \param L Length of the box
 */
/*template<typename T, typename U>
void pbcOffset(T &x, U &o, const T L) {
	o = (U) std::floor(x / L);
	x -= L * o;
}*/

#endif // ACTIVEBROWNIAN_STATE_H_
