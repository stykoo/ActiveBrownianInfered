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
#include <forward_list>
#include <random>

//! Type name for vector of positions
typedef std::vector< std::array<double, 2> > PositionVec;

class State; // Forward declaration

/*!
 * \brief Class for the boxes in which particles are.
 *
 * We divide the space into boxes of size approximately 1 and classify
 * the particles accordingly. This enables us to compute the internal forces
 * in (expected) linear time with respect to the number of particles.
 */
class Boxes {
	public:
		Boxes(const double len, const long n_parts);
		void update(const State *state); //!< Update according to state

		//!< Return the number of boxes
		long getNBoxes() const {
			return n_boxes;
		}

		//!< Return an iterator to the begin of the neighbors' list of box i.
		auto getNbrsBegin(const size_t i) const {
			return nbrs[i].cbegin();
		}
		//!< Return an iterator to the begin of the neighbors' list of box i.
		auto getNbrsEnd(const size_t i) const {
			return nbrs[i].cend();
		}

		//!< Return the box in which particle i is.
		long getBoxOfPart(const size_t i) const {
			return box_of_part[i];
		}

		//!< Return an iterator to the begin of the list of particles in box i.
		auto getPartsOfBoxBegin(const size_t i) const {
			return parts_of_box[i].cbegin();
		}
		//!< Return an iterator to the end of the list of particles in box i.
		auto getPartsOfBoxEnd(const size_t i) const {
			return parts_of_box[i].cend();
		}

	private:
		void computeNbrs(); //!< Compute the neighboring boxes of a given box

		const long n_boxes_x; //!< Number of boxes in one direction
		const long n_boxes; //!< Total number of boxes
		const double len_box; //!< Length of a box (approximately 1)
		const long n_parts; //!< Number of particles
		//! Neighboring boxes of a given box
		std::vector< std::array<long, 9> > nbrs; 

		std::vector<long> box_of_part; // Box in which each particle is
		//!< Particles in a given box
		std::vector< std::forward_list<long> > parts_of_box;
};

/*!
 * \brief Class for the state of the system
 *
 * This class takes care of the initialization
 * and the evolution of the state of the system.
 */
class State {
	public:
		//! Constructor of State
		State(const double _len, const long _n_parts,
		      const double _pot_strength, const double _temperature,
			  const double _rot_dif, const double _activity, const double _dt);
		void evolve(); //!< Do one time step

		//! Get the x coordinate of the position of particle i  
		double getPosX(size_t i) const {
			return positions[i][0];
		}
		//! Get the y coordinate of the position of particle i  
		double getPosY(size_t i) const {
			return positions[i][1];
		}

		//! Get angle of particle i
		double getAngle(size_t i) const {
			return angles[i];
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
		std::normal_distribution<double> noiseAngle;

		Boxes boxes; //!< Boxes for algorithm

		//! Positions of the particles
		std::vector< std::array<double, 2> > positions;
		std::vector<double> angles; //<! Angles
		std::vector< std::array<double, 2> > forces;  //!< Internal forces
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

/*! 
 * \brief Periodic boundary conditions on a segment (symmetric)
 * 
 * Update x to be between -L/2 and L/2.
 *
 * \param x Value
 * \param L Length of the box
 */
template<typename T>
void pbcSym(T &x, const T L) {
	x -= L * std::round(x / L);
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
template<typename T, typename U>
void pbcOffset(T &x, U &o, const T L) {
	o = (U) std::floor(x / L);
	x -= L * o;
}

/*! 
 * \brief Compute a to the power b for b positive integer
 * 
 * \param a Number
 * \param b Power (positive integer)
 * \return a to the power b
 */
template<typename T, typename U>
int mypow(const T a, const U b) {
    if (b <= 0)
        return 1;
    if (b % 2 == 1)
        return a * mypow(a, b-1);
    T c = mypow(a, b/2);
    return c * c;
}

#endif // ACTIVEBROWNIAN_STATE_H_
