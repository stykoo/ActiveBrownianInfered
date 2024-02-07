#include <string>
#include <exception>
#include <cassert>
#include "simul.h"
#include "state.h"

#ifndef NOVISU
#include <thread>
#include "visu/visu.h"
#endif

/* Create the simulation from parameter file. */
Simul::Simul(std::string fname) {
	status = SIMUL_INIT_SUCCESS;

	try {
		H5::H5File file(fname, H5F_ACC_RDONLY);
		loadAttribute(&file, &Lx, "Lx");
		loadAttribute(&file, &Ly, "Ly");
		loadAttribute(&file, &n_parts, "N");
		loadAttribute(&file, &diam, "d");
		loadAttribute(&file, &activity, "U");
		loadAttribute(&file, &trans_dif, "Dt");
		loadAttribute(&file, &rot_dif, "Dr");
		loadAttribute(&file, &pot_strength, "pot_strength");
		loadAttribute(&file, &dt, "dt");
		loadAttribute(&file, &n_iters, "niters");
		loadAttribute(&file, &n_iters_th, "niters_th");
		loadAttribute(&file, &skip, "skip");
		loadAttribute(&file, &dx, "dx");
		loadAttribute(&file, &r0, "r0");
		loadKs(file);
		loadCoeffs(file);
	} catch (H5::Exception& err) {
        err.printErrorStack();
		status = SIMUL_INIT_FAILED;
		return;
	}

	// Check if the values of the parameters are allowed
	if (
			notStrPos(Lx, "Lx")  || notStrPos(Ly, "Ly")
			|| notStrPos(n_parts, "n_parts")|| notStrPos(n_parts, "diam")
			|| notPos(activity, "U") || notPos(trans_dif, "Dt")
			|| notPos(rot_dif, "Dr") || notStrPos(dt, "dt")
			|| notPos(n_iters, "n_iters") || notPos(n_iters, "n_iters_th")
			|| notStrPos(skip, "skip") || notStrPos(dx, "dx")
			|| notStrPos(skip, "n_funs") || notStrPos(r0, "r0")
		) {
		status = SIMUL_INIT_FAILED;
		return;
	}

	// Structure for computation of infered interactions
	infered = new Infered(ks, coeffs, r0);
	
	// For correlations
	xmax = n_funs * r0;
	long n_div_pos = (long) xmax / dx;
	xmax = dx * n_div_pos;
	n_div = 2 * n_div_pos;
	n_frames_correl = 0;
	correls.assign(n_div * n_div, 0);

	/*for (int i = 0 ; i < 5 ; ++i) {
		trajectories.emplace_back();
		trajectories.back().assign(3 * n_parts, (double) i);
	}*/

#ifndef NOVISU
	sleep = 100;
#endif
}

Simul::~Simul() {
	delete infered;
}

void Simul::loadKs(H5::H5File &file) {
	H5::DataSet dset = file.openDataSet("ks");
	H5::DataSpace dspace = dset.getSpace();

	hsize_t rank;
	hsize_t dims[2];
	rank = dspace.getSimpleExtentDims(dims);
	assert(rank == 2 && dims[0] == 2);
	n_modes_tot = dims[1];

	H5::DataSpace memspace(1, &dims[1]);
	hsize_t dataCount[2] = {1, dims[1]};
	hsize_t dataOffset[2] = {0, 0};

	for (int i = 0 ; i < 2 ; ++i) {
		dataOffset[0] = i;
		ks[i].resize(n_modes_tot);
		dspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
		dset.read(ks[i].data(), H5::PredType::NATIVE_INT, memspace, dspace);
	}

	/*for (int i = 0 ; i < n_modes_tot ; ++i) {
		std::cout << ks[0][i] << " " << ks[1][i] << "\n";
	}
	std::cout << "\n";*/
}

void Simul::loadCoeffs(H5::H5File &file) {
	H5::DataSet dset = file.openDataSet("coeffs");
	H5::DataSpace dspace = dset.getSpace();

	hsize_t rank;
	hsize_t dims[3];
	rank = dspace.getSimpleExtentDims(dims);
	assert(rank == 3 && dims[0] == 3 && dims[1] == (hsize_t) n_modes_tot);

	n_funs = dims[2];

	hsize_t dimsm[2] = {dims[1], dims[2]};
	H5::DataSpace memspace(2, dimsm);
	hsize_t dataCount[3] = {1, dims[1], dims[2]};
	hsize_t dataOffset[3] = {0, 0, 0};

	for (int i = 0 ; i < 3 ; ++i) {
		dataOffset[0] = i;
		coeffs[i].resize(n_modes_tot * n_funs);
		dspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
		dset.read(coeffs[i].data(), H5::PredType::NATIVE_DOUBLE,
				  memspace, dspace);
	}

	/*for (int i = 0 ; i < 3 ; ++i) {
		for (int j = 0 ; j < n_modes_tot ; ++j) {
			for (int k = 0 ; k < n_funs ; ++k) {
				std::cout << coeffs[i][j*n_funs+k] << " ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}*/
}

/* Print the parameters of the simulation. */
void Simul::print() const {
	std::cout << "# [ActiveBrownianInfered] ";
	printParam(Lx, "Lx");
	printParam(Ly, "Ly");
	printParam(n_parts, "N");
	printParam(diam, "d");
	printParam(activity, "U");
	printParam(trans_dif, "Dt");
	printParam(rot_dif, "Dr");
	printParam(dt, "dt");
	printParam(n_iters, "niters");
	printParam(n_iters_th, "niters_th");
	printParam(skip, "skip");
	printParam(dx, "dx");
	printParam(n_funs, "nfuns");
	printParam(r0, "r0");
	std::cout << std::endl;
}

/* Run the simulation. */
void Simul::run() {
	if (status != SIMUL_INIT_SUCCESS) {
		std::cerr << "You should not be runing a failed simulation..."
		          << std::endl;
		return;
	}
	
	// Initialize the state of the system
	State state(Lx, Ly, n_parts, diam, trans_dif, rot_dif, activity,
				dt, pot_strength, *infered);
	
#ifndef NOVISU
	// Start thread for visualization
	Visu visu(&state, Lx, Ly, n_parts, diam / 2.);
	std::thread thVisu(&Visu::run, &visu); 
#endif

	// Thermalization
	for (long t = 0 ; t < n_iters_th ; ++t) {
		state.evolve();
	}
	// Time evolution
	for (long t = 0 ; t < n_iters ; ++t) {
#ifndef NOVISU
		std::cout << " t = " << t * dt << "\r" << std::flush;
#endif
		if (t % skip == 0) {
			// Compute correlations
			addCorrelations(&state);
			// Save trajectories
			state.store(trajectories.emplace_back());
		}
		state.evolve();
#ifndef NOVISU
		std::this_thread::sleep_for(std::chrono::milliseconds(SLEEP));
#endif
	}

#ifndef NOVISU
		thVisu.join();
#endif
}

/* Add data to the correlations (at current time) */
void Simul::addCorrelations(const State *state) {
	n_frames_correl++;

	const auto& pos_x = state->getPosX();
	const auto& pos_y = state->getPosY();
	const auto& angles = state->getAngles();

	// Non-MKL version
	// For each pair of particles
	for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = i + 1 ; j < n_parts ; ++j) {
			double x = pos_x[j] - pos_x[i];
			double y = pos_y[j] - pos_y[i];
			pbcSym(x, Lx);
			pbcSym(y, Ly);
			double r = std::sqrt(x * x + y * y);
			if (r >= xmax) { // Get rid of the points too far away
				continue;
			}

			double phi = std::atan2(y, x);
			double theta1 = angles[j] - phi;

			double xr = std::cos(theta1) * r;
			double yr = std::sin(theta1) * r;
			size_t b1 = (size_t) ((xmax + xr) / dx);
			size_t b2 = (size_t) ((xmax + yr) / dx);
			correls[b1 * n_div + b2]++; // Add 1 in the right box
		}
	}

/*
	// MKL version
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
		size_t b1 = (size_t) dxs[k];
		size_t b2 = (size_t) dys[k];
		correls[b1 * n_div_r + b2]++; // Add 1 in the right box
	}
*/
}

/* Save the simulation */
void Simul::save(std::string fname) {
	std::cout << "Saving simulation to " << fname << "\n";
	try {
		H5::H5File file(fname, H5F_ACC_TRUNC);

		writeAttribute(&file, H5::PredType::NATIVE_DOUBLE, &Lx, "Lx");
		writeAttribute(&file, H5::PredType::NATIVE_DOUBLE, &Ly, "Ly");
		writeAttribute(&file, H5::PredType::NATIVE_LONG, &n_parts, "N");
		writeAttribute(&file, H5::PredType::NATIVE_DOUBLE, &diam, "d");
		writeAttribute(&file, H5::PredType::NATIVE_DOUBLE, &activity, "U");
		writeAttribute(&file, H5::PredType::NATIVE_DOUBLE, &trans_dif, "Dt");
		writeAttribute(&file, H5::PredType::NATIVE_DOUBLE, &rot_dif, "Dr");
		writeAttribute(&file, H5::PredType::NATIVE_DOUBLE, &dt, "dt");
		writeAttribute(&file, H5::PredType::NATIVE_LONG, &n_iters, "niters");
		writeAttribute(&file, H5::PredType::NATIVE_LONG, &n_iters_th,
				       "niters_th");
		writeAttribute(&file, H5::PredType::NATIVE_LONG, &skip, "skip");
		writeAttribute(&file, H5::PredType::NATIVE_DOUBLE, &dx, "dx");
		writeAttribute(&file, H5::PredType::NATIVE_DOUBLE, &r0, "r0");

		writeKs(file);
		writeCoeffs(file);
		computeAndWriteForces(file);
		writeTrajectories(file);
		writeCorrelations(file);
	} catch (H5::Exception& err) {
        err.printErrorStack();
	}
}

void Simul::writeKs(H5::H5File &file) {
	hsize_t dims[2] = {2, (hsize_t) n_modes_tot};
	H5::DataSpace dspace(2, dims);
	H5::DataSet dset = file.createDataSet("ks", H5::PredType::NATIVE_INT,
							              dspace);

	H5::DataSpace memspace(1, &dims[1]);
	hsize_t dataCount[2] = {1, dims[1]};
	hsize_t dataOffset[2] = {0, 0};

	for (int i = 0 ; i < 2 ; ++i) {
		dataOffset[0] = i;
		dspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
		dset.write(ks[i].data(), H5::PredType::NATIVE_INT, memspace, dspace);
	}
}

void Simul::writeCoeffs(H5::H5File &file) {
	hsize_t dims[3] = {3, (hsize_t) n_modes_tot, (hsize_t) n_funs};
	H5::DataSpace dspace(3, dims);
	H5::DataSet dset = file.createDataSet("coeffs",
			                              H5::PredType::NATIVE_DOUBLE, dspace);

	hsize_t dimsm[2] = {dims[1], dims[2]};
	H5::DataSpace memspace(2, dimsm);
	hsize_t dataCount[3] = {1, dims[1], dims[2]};
	hsize_t dataOffset[3] = {0, 0, 0};

	for (int i = 0 ; i < 3 ; ++i) {
		dataOffset[0] = i;
		dspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
		dset.write(coeffs[i].data(), H5::PredType::NATIVE_DOUBLE, memspace,
				   dspace);
	}
}

void Simul::computeAndWriteForces(H5::H5File &file) {
	double rmax = n_funs * r0;
	long n_div_r = (long) rmax / dx;
	double step_angle = (2 * M_PI) / N_DIV_ANGLE;

	// New group
	H5::Group group = file.createGroup("/force_fields");

	// Compute
	std::vector<double> rlist(n_div_r), thlist(N_DIV_ANGLE);
	std::array<std::vector<double>, 3> forces;
	for (int i = 0 ; i < 3 ; ++i) {
		forces[i].resize(N_DIV_ANGLE*N_DIV_ANGLE*n_div_r);
	}

	for (long i = 0 ; i < n_div_r ; ++i) {
		rlist[i] = i * dx;
	}
	for (long j = 0 ; j < N_DIV_ANGLE ; ++j) {
		thlist[j] = j * step_angle;
	}

	// Write
	hsize_t dims_r[1] = {(hsize_t) n_div_r,};
	H5::DataSpace dspace_r(1, dims_r);
	H5::DataSet dset_r = file.createDataSet("force_fields/r_list",
			                                H5::PredType::NATIVE_DOUBLE,
						  				    dspace_r);
	dset_r.write(rlist.data(), H5::PredType::NATIVE_DOUBLE);

	hsize_t dims_t[1] = {N_DIV_ANGLE, };
	H5::DataSpace dspace_t(1, dims_t);
	H5::DataSet dset_t = file.createDataSet("force_fields/th_list",
			                                H5::PredType::NATIVE_DOUBLE,
						  				    dspace_t);
	dset_t.write(thlist.data(), H5::PredType::NATIVE_DOUBLE);

	// Compute
	double fr, ft, fo;
	long b;
	for (long i = 0 ; i < n_div_r ; ++i) {
		for (long j = 0 ; j < N_DIV_ANGLE ; ++j) {
			for (long k = 0 ; k < N_DIV_ANGLE ; ++k) {
				infered->computeForces(rlist[i], thlist[j], thlist[k],
						               fr, ft, fo);
				b = i*N_DIV_ANGLE*N_DIV_ANGLE + j*N_DIV_ANGLE + k;
				forces[0][b] = fr;
				forces[1][b] = ft;
				forces[2][b] = fo;
			}
		}
	}

	// Write
	hsize_t dims[4] = {3, (hsize_t) n_div_r, (hsize_t) N_DIV_ANGLE,
		               (hsize_t) N_DIV_ANGLE};
	H5::DataSpace dspace(4, dims);
	H5::DataSet dset = file.createDataSet("force_fields/forces",
			                               H5::PredType::NATIVE_DOUBLE,
										   dspace);

	hsize_t dimsm[3] = {dims[1], dims[2], dims[3]};
	H5::DataSpace memspace(3, dimsm);
	hsize_t dataCount[4] = {1, dims[1], dims[2], dims[3]};
	hsize_t dataOffset[4] = {0, 0, 0, 0};

	for (int i = 0 ; i < 3 ; ++i) {
		dataOffset[0] = i;
		dspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
		dset.write(forces[i].data(), H5::PredType::NATIVE_DOUBLE, memspace,
				   dspace);
	}
}

void Simul::writeTrajectories(H5::H5File &file) {
	H5::Group group = file.createGroup("/trajectories");

	double dt_traj = skip * dt;
	writeAttribute(&group, H5::PredType::NATIVE_DOUBLE, &dt_traj,
			       "dt");

	// Correlations
	size_t n_frames = trajectories.size();
	hsize_t dims[3] = {n_frames, 3, (hsize_t) n_parts};
	H5::DataSpace dspace(3, dims);
	
	hsize_t chunk_dims[3] = {1, 1, (hsize_t) n_parts};
	H5::DSetCreatPropList plist;
	plist.setDeflate(6);
	plist.setChunk(3, chunk_dims);

	H5::DataSet dset = \
		file.createDataSet("trajectories/traj", H5::PredType::NATIVE_DOUBLE,
						   dspace, plist);

	H5::DataSpace memspace(2, &dims[1]);
	hsize_t dataCount[3] = {1, dims[1], dims[2]};
	hsize_t dataOffset[3] = {0, 0, 0};

	for (size_t i = 0 ; i < n_frames ; ++i) {
		dataOffset[0] = i;
		dspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
		dset.write(trajectories[i].data(), H5::PredType::NATIVE_DOUBLE,
				   memspace, dspace);
	}
}

void Simul::writeCorrelations(H5::H5File &file) {
	// New group
	H5::Group group = file.createGroup("/correlations");

	writeAttribute(&group, H5::PredType::NATIVE_LONG, &n_frames_correl,
			       "n_frames");

	// Edges of boxes
	std::vector<double> xedges(n_div+1);
	for (long i = 0 ; i <= n_div ; ++i) {
		xedges[i] = -xmax + i * dx;
	}

	hsize_t dims_x[1] = {(hsize_t) n_div,};
	H5::DataSpace dspace_x(1, dims_x);
	H5::DataSet dset_r = file.createDataSet("correlations/x_edges",
			                                H5::PredType::NATIVE_DOUBLE,
						  				    dspace_x);
	dset_r.write(xedges.data(), H5::PredType::NATIVE_DOUBLE);

	// Correlations
	hsize_t dims[2] = {(hsize_t) n_div, (hsize_t) n_div};
	H5::DataSpace dspace(2, dims);
	
	hsize_t chunk_dims[2] = {1, (hsize_t) std::min(1000l, n_div)};
	H5::DSetCreatPropList plist;
	plist.setDeflate(6);
	plist.setChunk(2, chunk_dims);

	H5::DataSet dset = \
		file.createDataSet("correlations/correl", H5::PredType::NATIVE_LLONG,
						   dspace, plist);
	dset.write(correls.data(), H5::PredType::NATIVE_LLONG);
}
