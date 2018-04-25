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
 * \file boxes.h
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Boxes for efficient algorithm.
 *
 * It defines and implement the class Boxes in arbitrary dimension.
 */

#ifndef ACTIVEBROWNIAN_BOXES_H_
#define ACTIVEBROWNIAN_BOXES_H_

#include <vector>
#include <array>

/*!
 * \brief Class for the boxes in which particles are.
 *
 * We divide the space into boxes of size approximately 1 and classify
 * the particles accordingly. This enables us to compute the internal forces
 * in (expected) linear time with respect to the number of particles.
 */
template<int DIM>
class Boxes {
	public:
		Boxes(const double len, const long n_parts, const int fac=1);
		//! Update according to the positions
		void update(const std::array< std::vector<double>, DIM> &pos);

		//! Return the number of boxes
		long getNBoxes() const {
			return n_boxes;
		}

		const std::vector< std::vector<long> > & getNbrsPos() const {
			return nbrs_pos;
		}

		const std::vector< std::vector<long> > & getPartsOfBox() const {
			return parts_of_box;
		}

	private:
		//!< Compute the neighboring boxes of a given box
		void computeNbrsPos();

		const long n_boxes_x; //!< Number of boxes in one direction
		const long n_boxes; //!< Total number of boxes
		const double len_box; //!< Length of a box (approximately 1/fac)
		const long n_parts; //!< Number of particles
		const int fac; //!< Factor for the size of the boxes
		//! Neighboring boxes of a given box along the 'positive' directions
		std::vector< std::vector<long> > nbrs_pos; 

		//!< Particles in a given box
		std::vector< std::vector<long> > parts_of_box;
};

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

/*
 * \brief Constructor of boxes.
 *
 * Construct the boxes with which we will be able to classify the particles.
 *
 * \param len Length of the system
 * \param n_parts Number of particles
 */
template<int DIM>
Boxes<DIM>::Boxes(const double len, const long n_parts, const int fac) :
		n_boxes_x(fac * std::floor(len)), n_boxes(mypow(n_boxes_x, DIM)),
		len_box(len / n_boxes_x), n_parts(n_parts), fac(fac) {
	nbrs_pos.resize(n_boxes);
	computeNbrsPos();
	parts_of_box.resize(n_boxes);
}

/*
 * \brief Classify the particles at given positions in the boxes.
 * 
 * Warning: the positions should be consistent with the Boxes object!
 *
 * \param pos Positions of the particles
 */
template<int DIM>
void Boxes<DIM>::update(const std::array<std::vector<double>, DIM> &pos) {
    for (long i=0 ; i < n_boxes ; ++i) {
		parts_of_box[i].clear();
	}

    for (long i=0 ; i < n_parts ; ++i) {
        long box = 0;
        for (int a = 0 ; a < DIM ; a++) {
			//long ba = (long) std::floor((*pos)[i][a] / len_box);
			// Round towards 0
			long ba = (long) (pos[a][i] / len_box);
            box += mypow(n_boxes_x, a) * ba;
        }

		parts_of_box[box].push_back(i);
	}
}

/*
 * \brief Classify the particles at given positions in the boxes.
 * 
 * Specialization for d = 2 (removes the call to mypow).
 * Warning: the positions should be consistent with the Boxes object!
 *
 * \param pos Positions of the particles
 */
template<>
void Boxes<2>::update(const std::array<std::vector<double>, 2> &pos) {
    for (long i=0 ; i < n_boxes ; ++i) {
		parts_of_box[i].clear();
	}

    for (long i=0 ; i < n_parts ; ++i) {
        long box = 0;
		// Round towards 0
		long bx = (long) (pos[0][i] / len_box);
		long by = (long) (pos[1][i] / len_box);
		box = bx + n_boxes_x * by;

		parts_of_box[box].push_back(i);
	}
}

/*
 * \brief Compute the indices of the neighboring boxes of each box
 * along the 'positive' direction of each axis.
 *
 * This is done only once and does not need to be efficient.
 */
template<int DIM>
void Boxes<DIM>::computeNbrsPos() {
	std::array<double, DIM> pows_n_boxes_x;
	for (int a = 0 ; a < DIM ; ++a) {
		pows_n_boxes_x[a] = mypow(n_boxes_x, a);
	}
	std::array<long, DIM> coos;

	// for all the particles
	for (long k = 0 ; k < n_boxes ; ++k) {
		// Coordinates of the box
		long i = k;
		for (int a = DIM-1 ; a >= 0 ; --a) {
			coos[a] = i / pows_n_boxes_x[a];
			i -= pows_n_boxes_x[a] * coos[a];
		}

		nbrs_pos[k].clear();
		nbrs_pos[k].push_back(k);

		// First dimension (half space)
		for (int b = 0 ; b < fac ; ++b) {
			long d = (coos[0] + 1 + b) % n_boxes_x - coos[0];
			nbrs_pos[k].push_back(k + d);

			// Second dimension for k + d (both sides)
			for (int c = 0 ; c < fac ; ++c) {
				long e1 = (coos[1] + c + 1) % n_boxes_x - coos[1];
				long e2 = (coos[1] - c - 1 + n_boxes_x) % n_boxes_x - coos[1];
				nbrs_pos[k].push_back(k + d + n_boxes_x * e1);
				nbrs_pos[k].push_back(k + d + n_boxes_x * e2);
			}
		}
		// Second dimension for k (only half space)
		for (int b = 0 ; b < fac ; ++b) {
			long d = (coos[1] + b + 1) % n_boxes_x - coos[1];
			nbrs_pos[k].push_back(k + n_boxes_x * d);
		}

		// Other dimensions (both sides)
		for (int a = 2 ; a < DIM ; ++a) {
			std::vector<long> tmp = nbrs_pos[k];

			for (int c = 0 ; c < fac ; ++c) {
				long e1 = (coos[a] + c + 1) % n_boxes_x - coos[a];
				long e2 = (coos[a] - c - 1 + n_boxes_x) % n_boxes_x - coos[a];
				for (long x : tmp) {
					nbrs_pos[k].push_back(x + pows_n_boxes_x[a] * e1);
					nbrs_pos[k].push_back(x + pows_n_boxes_x[a] * e2);
				}
			}
		}
	}
}

#endif // ACTIVEBROWNIAN_BOXES_H_
