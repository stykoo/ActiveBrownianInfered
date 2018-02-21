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
 * \file visu.h
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Represent a point on an unit sphere.
 */

#ifndef ACTIVEBROWNIAN_POINTONSPHERE_H_
#define ACTIVEBROWNIAN_POINTONSPHERE_H_

#include <cmath>
#include <array>
#include <random>

/*!
 * \brief Class for a random point on a unit sphere.
 */
class PointOnSphere {
	public:
		PointOnSphere(std::mt19937 &rng);
		PointOnSphere(const PointOnSphere &p) {
			coos = p.coos;
		}

		void randomRotation(const double stddev, std::mt19937 &rng);
		void renormalize();
		
		double getX() const {
			return coos[0];
		}
		double getY() const {
			return coos[1];
		}
		double getZ() const {
			return coos[2];
		}
		double getNorm() const {
			// This should always be 1
			return std::sqrt(coos[0]*coos[0] + coos[1]*coos[1]
					         + coos[2]*coos[2]);
		}

		//! Return theta of spherical coordinates
		double getTheta() const {
			return std::acos(coos[2]);
		}
		//! Return phi of spherical coordinates
		double getPhi() const {
			if (coos[0] == 0.0 && coos[1] == 0.0) {
				return 0; // Arbitrary
			} else {
				return std::atan2(coos[1], coos[0]);
			}
		}

	private:
		std::array<double, 3> coos; //!< Coordinates
};

#endif // ACTIVEBROWNIAN_POINTONSPHERE_H_
