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
 * \file observables.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Computation and export of the observables
*/

#include <iostream>
//#include <cassert>
#include <cmath>
#include "H5Cpp.h"
#include "observables.h"

/*
 * \brief Constructor of Observables
 *
 * Initialize the vector for correlations.
 */
Observables::Observables(const double len_, const long n_parts_,
		                 const double step_r_) :
		len(len_), n_parts(n_parts_), step_r(step_r_), scal_r(1.0 / step_r)
#ifdef USE_MKL
		, n_pairs(n_parts * (n_parts - 1) / 2),
		dxs(n_pairs), dys(n_pairs),
		phis(n_pairs), drs(n_pairs), thetas1(n_pairs)
#endif
{
	n_div_r = (long) std::ceil(len / step_r);
	n_div_tot = n_div_r * n_div_r;
	step_r = len / n_div_r; // The step divides exactly the box
	scal_r = 1.0 / step_r;
	n_calls = 0;
	correls.assign(n_div_tot, 0);
}


/*
 * \brief Compute the observables for a given state
 */
void Observables::compute(const State *state) {
	const std::vector<double> & pos_x = state->getPosX();
	const std::vector<double> & pos_y = state->getPosY();
	const std::vector<double> & angles = state->getAngles();

	n_calls++;

#ifdef USE_MKL // MKL version
	long k = 0;
	for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = i + 1 ; j < n_parts ; ++j) {
			dxs[k] = pos_x[j] - pos_x[i];
			dys[k] = pos_y[j] - pos_y[i];
			thetas1[k] = angles[j];
			++k;
		}
	}

	pbcSymMKL(dxs, len, phis, n_pairs);
	pbcSymMKL(dys, len, phis, n_pairs);
	// Computation of phis
	vdAtan2(n_pairs, dys.data(), dxs.data(), phis.data());
	// Computation of drs
	vdSqr(n_pairs, dxs.data(), dxs.data());
	vdSqr(n_pairs, dys.data(), dys.data());
	vdAdd(n_pairs, dxs.data(), dys.data(), drs.data());
	vdSqrt(n_pairs, drs.data(), drs.data());
	cblas_dscal(n_pairs, scal_r, drs.data(), 1); // Scaling
	// Computation of thetas1
	cblas_daxpy(n_pairs, -1.0, phis.data(), 1, thetas1.data(), 1);
	// r cos(theta1), r sin(theta1)
	vdSinCos(n_pairs, thetas1.data(), dys.data(), dxs.data());
	vdMul(n_pairs, dxs.data(), drs.data(), dxs.data());
	vdMul(n_pairs, dys.data(), drs.data(), dys.data());
	pbcMKL(dxs, n_div_r, phis, n_pairs);
	pbcMKL(dys, n_div_r, phis, n_pairs);

	double rmax = len * scal_r / 2;
	for (long k = 0 ; k < n_pairs ; ++k) {
		if (drs[k] >= rmax) { // Get rid of the points too far away
			continue;
		}
		size_t box, b1, b2;
		b1 = (size_t) dxs[k];
		b2 = (size_t) dys[k];
		box = b1 * n_div_r + b2;
		correls[box]++; // Add 1 in the right box
	}
#else // Basic version
	// For each pair of particles
	for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = i + 1 ; j < n_parts ; ++j) {
			double dx = pos_x[j] - pos_x[i];
			pbcSym(dx, len);
			double dy = pos_y[j] - pos_y[i];
			pbcSym(dy, len);
			double dr = std::sqrt(dx * dx + dy * dy);
			if (dr > len / 2.) { // Get rid of the points too far away
				continue;
			}

			double phi = std::atan2(dy, dx);
			double theta1 = angles[j] - phi;
			pbc(theta1, 2 * M_PI);

			size_t box = 0;
			double x = std::cos(theta1) * dr;
			double y = std::sin(theta1) * dr;
			pbc(x, len);
			pbc(y, len);
			size_t b1 = (size_t) (x * scal_r);
			size_t b2 = (size_t) (y * scal_r);
			box = b1 * n_div_r + b2;

			correls[box]++; // Add 1 in the right box
		}
	}
#endif
}

/*
 * \brief Export the observables to a hdf5 file
 */
void Observables::writeH5(const std::string fname, double rho, long n_parts,
				          double temperature, double rot_dif,
						  double activity, double dt,
						  long n_iters, long n_iters_th, long skip) const {
	try {
		H5::H5File file(fname, H5F_ACC_TRUNC);

		// General attributes
		H5::DataSpace default_ds;
		H5::Attribute a_rho = file.createAttribute(
				"rho", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_rho.write(H5::PredType::NATIVE_DOUBLE, &rho);
		H5::Attribute a_n_parts = file.createAttribute(
				"n_parts", H5::PredType::NATIVE_LONG, default_ds);
		a_n_parts.write(H5::PredType::NATIVE_LONG, &n_parts);
		H5::Attribute a_temperature = file.createAttribute(
				"temperature", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_temperature.write(H5::PredType::NATIVE_DOUBLE, &temperature);
		H5::Attribute a_rot_dif = file.createAttribute(
				"rot_dif", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_rot_dif.write(H5::PredType::NATIVE_DOUBLE, &rot_dif);
		H5::Attribute a_activity = file.createAttribute(
				"activity", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_activity.write(H5::PredType::NATIVE_DOUBLE, &activity);
		H5::Attribute a_dt = file.createAttribute(
				"dt", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_dt.write(H5::PredType::NATIVE_DOUBLE, &dt);
		H5::Attribute a_n_iters = file.createAttribute(
				"n_iters", H5::PredType::NATIVE_LONG, default_ds);
		a_n_iters.write(H5::PredType::NATIVE_LONG, &n_iters);
		H5::Attribute a_n_iters_th = file.createAttribute(
				"n_iters_th", H5::PredType::NATIVE_LONG, default_ds);
		a_n_iters_th.write(H5::PredType::NATIVE_LONG, &n_iters_th);
		H5::Attribute a_skip = file.createAttribute(
				"skip", H5::PredType::NATIVE_LONG, default_ds);
		a_skip.write(H5::PredType::NATIVE_LONG, &skip);
		
		// We chunk the data and compress it
		// Chunking should depend on how we intend to read the data
		hsize_t chunk_dims[3]; // Dimension of a chunk (compression)
		hsize_t dims[3]; // Dimensions of the data
		int ndim = 3;
		ndim = 2;
		chunk_dims[0] = 1;
		chunk_dims[1] = std::min(1000l, n_div_r);
		dims[0] = (hsize_t) n_div_r;
		dims[1] = (hsize_t) n_div_r;

		// Create dataset
		H5::DSetCreatPropList plist;
		plist.setDeflate(6);
		plist.setChunk(ndim, chunk_dims);
		H5::DataSpace dataspace(ndim, dims);
		H5::DataSet dataset(
			file.createDataSet("correlations", H5::PredType::NATIVE_LLONG,
							   dataspace, plist)
			);
		// Write data
		dataset.write(correls.data(), H5::PredType::NATIVE_LLONG);

		// Attributes for correlations
		H5::Attribute a_dr = dataset.createAttribute(
				"dr", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_dr.write(H5::PredType::NATIVE_DOUBLE, &step_r);
	} catch (H5::Exception& err) {
        err.printErrorStack();
	}
}
