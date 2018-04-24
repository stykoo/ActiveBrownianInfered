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
		Boxes(const double len, const long n_parts);
		//! Update according to the positions
		void update(const std::vector< std::array<double, DIM> > *pos);

		//! Return the number of boxes
		long getNBoxes() const {
			return n_boxes;
		}

		const std::vector< std::vector<long> > * getNbrsPos() const {
			return &nbrs_pos;
		}

		const std::vector< std::vector<long> > * getPartsOfBox() const {
			return &parts_of_box;
		}

	private:
		//!< Compute the neighboring boxes of a given box
		void computeNbrsPos();

		const long n_boxes_x; //!< Number of boxes in one direction
		const long n_boxes; //!< Total number of boxes
		const double len_box; //!< Length of a box (approximately 1)
		const long n_parts; //!< Number of particles
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
Boxes<DIM>::Boxes(const double len, const long n_parts) :
		n_boxes_x(std::floor(len)), n_boxes(mypow(n_boxes_x, DIM)),
		len_box(len / n_boxes_x), n_parts(n_parts) {
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
void Boxes<DIM>::update(const std::vector< std::array<double, DIM> > *pos) {
    for (long i=0 ; i < n_boxes ; ++i) {
		parts_of_box[i].clear();
	}

    for (long i=0 ; i < n_parts ; ++i) {
        long box = 0;
        for (int a = 0 ; a < DIM ; a++) {
			//long ba = (long) std::floor((*pos)[i][a] / len_box);
			// Round towards 0
			long ba = (long) ((*pos)[i][a] / len_box);
            box += mypow(n_boxes_x, a) * ba;
        }

		// box_of_part[i] = box;
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
void Boxes<2>::update(const std::vector< std::array<double, 2> > *pos) {
    for (long i=0 ; i < n_boxes ; ++i) {
		parts_of_box[i].clear();
	}

    for (long i=0 ; i < n_parts ; ++i) {
        long box = 0;
		// Round towards 0
		long bx = (long) ((*pos)[i][0] / len_box);
		long by = (long) ((*pos)[i][1] / len_box);
		box = bx + n_boxes_x * by;

		//box_of_part[i] = box;
		parts_of_box[box].push_back(i);
	}
}

/*
 * \brief Compute the indices of the neighboring boxes of each box
 * along the 'positive' direction of each axis.
 */
template<int DIM>
void Boxes<DIM>::computeNbrsPos() {
	for (long k = 0 ; k < n_boxes ; ++k) {
		nbrs_pos[k].clear();
		nbrs_pos[k].push_back(k);

		long i = k;

		for (int a = DIM-1 ; a >= 0 ; a--) {
			std::vector<long> tmp = nbrs_pos[k];

			long p = mypow(n_boxes_x, a);
			long ia = i / p; // Coordinate along axis a
			// Periodic shift of +1 along axis a
			long d = (ia + 1) % n_boxes_x - ia;
			i -= p * ia;

			// We add shift to the points we have so far
			for (long x : tmp) {
				nbrs_pos[k].push_back(x + p * d);
			}
		}
	}
}

#endif // ACTIVEBROWNIAN_BOXES_H_
