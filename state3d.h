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
 * \brief State of the system in dimension 3
 *
 * Header file for state3d.cpp
 * It defines the class State3d.
 */

#ifndef ACTIVEBROWNIAN_STATE3D_H_
#define ACTIVEBROWNIAN_STATE3D_H_

#include <vector>
#include <array>
#include <random>
#include "pointOnSphere.h"
#include "boxes.h"

/*!
 * \brief Class for the state of the system in dimension 3
 *
 * This class takes care of the initialization
 * and the evolution of the state of the system.
 */
class State3d {
	public:
		//! Constructor of State
		State3d(const double _len, const long _n_parts,
		        const double _pot_strength, const double _temperature,
			    const double _rot_dif, const double _activity,
				const double _dt);
		void evolve(); //!< Do one time step

		//! Get the x coordinate of the position of particle i  
		double getPosX(size_t i) const {
			return positions[i][0];
		}
		//! Get the y coordinate of the position of particle i  
		double getPosY(size_t i) const {
			return positions[i][1];
		}
		//! Get the z coordinate of the position of particle i  
		double getPosZ(size_t i) const {
			return positions[i][2];
		}

		//! Get orientation of particle i
		PointOnSphere getOrient(size_t i) const {
			return orients[i];
		}


	private:
		void calcInternalForces(); //!< Compute internal forces
		void enforcePBC(); //!< Enforce periodic boundary conditions

		const double len; //!< Length of the box
		const long n_parts; //!< Number of particles
		const double pot_strength; //!< Strength of the interparticle potential
		const double activity; //!< Activity
		const double dt; //!< Timestep

		std::mt19937 rng; //!< Random number generator
		//! Gaussian noise for temperature
		std::normal_distribution<double> noiseTemp;
		//! Gaussian noise for angle
		const double stddevOrient;

		Boxes<3> boxes; //!< Boxes for algorithm

		//! Positions of the particles
		std::vector< std::array<double, 3> > positions;
		std::vector<PointOnSphere> orients; //<! Orientations
		std::vector< std::array<double, 3> > forces;  //!< Internal forces
};

#endif // ACTIVEBROWNIAN_STATE3D_H_
