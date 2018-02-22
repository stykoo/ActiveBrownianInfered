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
 * \file simul.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Represent a point on an unit sphere.
 */

#include "pointOnSphere.h"

PointOnSphere::PointOnSphere(std::mt19937 &rng) {
	// Random point on a sphere
	// From W. Krauth's book
	std::normal_distribution<double> distrib(0.0, 1.0);
	double x = distrib(rng);
	double y = distrib(rng);
	double z = distrib(rng);
	double n = std::sqrt(x * x + y * y + z * z);
	coos[0] = x / n;
	coos[1] = y / n;
	coos[2] = z / n;
}

void PointOnSphere::randomRotation(const double stddev, std::mt19937 &rng) {
	// Find u and v to complete an orthonormal basis
	// See https://graphics.pixar.com/library/OrthonormalB/paper.pdf
	std::array<double, 3> u, v;
	double s = std::copysign(1.0, coos[2]);
	double a = -1.0 / (s + coos[2]);
	double b = coos[0] * coos[1] * a;
	u[0] = 1.0 + s * coos[0] * coos[0] * a;
	u[1] = s * b;
	u[2] = -s * coos[0];
	v[0] = b;
	v[1] = s * coos[1] * coos[1] * a;
	v[2] = -coos[1];

	// Draw a uniform angle and a gaussian angle
	std::uniform_real_distribution<double> distUnif(0, 2.0 * M_PI);
	std::normal_distribution<double> distGauss(0, stddev);
	double phi = distUnif(rng);
	double theta = distGauss(rng);

	// Generate new point
	double cr = std::cos(theta);
	double st = std::sin(theta);
	double cu = st * std::cos(phi);
	double cv = st * std::sin(phi);
	coos[0] = cr*coos[0] + cu*u[0] + cv*u[0];
	coos[1] = cr*coos[1] + cu*u[1] + cv*u[1];
	coos[2] = cr*coos[2] + cu*u[2] + cv*u[2];
}

// In priciple this should never be needed.
void PointOnSphere::renormalize() {
	double n = getNorm();
	coos[0] /= n;
	coos[1] /= n;
	coos[2] /= n;
}
