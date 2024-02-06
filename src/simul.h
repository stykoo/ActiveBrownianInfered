#ifndef ACTIVEBROWNIAN_SIMUL_H_
#define ACTIVEBROWNIAN_SIMUL_H_

#include <array>
#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include "H5Cpp.h"
#include "infered.h"
#include "state.h"

#define SLEEP 100 // Number of milliseconds before updating visualization
#define N_DIV_ANGLE 50 // Number of angular divisions

//! State of the simulation after initialization
enum SimulInitStatus {
	SIMUL_INIT_SUCCESS, //!< Successful initialization
	SIMUL_INIT_FAILED //!< Failed initialization
};

/* Class for simulation.  */
class Simul {
	public:
		Simul(std::string fname); //!< Constructor from file name
		~Simul();

		void print() const; //!< Print the parameters
		void run(); //!< Run the simulation
		void save(std::string fname); //!< Save the simulation

		//! Get initialization status
		SimulInitStatus getStatus() const { return status; }

	private:
		void loadKs(H5::H5File &file);
		void loadCoeffs(H5::H5File &file);
		void addCorrelations(const State *state);
		void writeKs(H5::H5File &file);
		void writeCoeffs(H5::H5File &file);
		void computeAndWriteForces(H5::H5File &file);
		void writeTrajectories(H5::H5File &file);
		void writeCorrelations(H5::H5File &file);

		// NEW
		SimulInitStatus status; //!< Status after initialization
								
		double Lx, Ly; //!< Size of the box
		long n_parts; //!< Number of particles
		double diam; //!< Particle diameter
		double activity; //!< Activity
		double trans_dif; //!< Translational diffusivity
		double rot_dif; //!< Rotational diffusivity
		double dt; //!< Timestep
		double pot_strength; //!< Stength of repulsive potential
		long n_iters; //!< Number of time iterations
		long n_iters_th; //!< Number of time iterations of thermalization
		long skip; //!< Iterations between two computation of observables
		double dx; //!< Spatial resolution for correlations
		double r0; //!< Radial increment for radial functions
		int n_funs; //!< Number of radial functions
		int n_modes_tot; //! Total number of angular modes
		std::array<std::vector<int>, 2> ks; //!< Factors for cos / sin
		std::array<std::vector<double>, 3> coeffs; //!< Infered coefficients
		
		Infered *infered; // Inference

		std::deque<std::vector<double>> trajectories;

		// Stuff for correlations
		double xmax;
		long n_div;
		long n_frames_correl;
		std::vector<long long> correls; //!< Correlations

#ifndef NOVISU
		int sleep; //!< Number of milliseconds to sleep for between iterations
#endif
};

/* Load attribute from HDF5 file or group. */
template<typename T>
void loadAttribute(H5::Group *g, T *data, std::string name) {
	H5::Attribute a = g->openAttribute(name);
	a.read(a.getDataType(), data); 
}

/* Write attribute to HDF5 file or group. */
template<typename T>
void writeAttribute(H5::Group *g, const H5::DataType &type, T *data,
		            std::string name) {
	H5::DataSpace dspace;
	H5::Attribute a = g->createAttribute(name, type, dspace);
	a.write(type, data);
}

/* Print a parameter. */
template<typename T> void printParam(const T &a, std::string name) {
	std::cout << name << "=" << a << ", ";
}

/* Check if variable is positive or zero. */
template<typename T> bool notPos(const T &a, std::string name) {
	if (a < T(0)) {
		std::cerr << "Error: " << name << " should be positive."
		          << std::endl;
		return true;
	}
	return false;
}

/* Check if variable is strictly positive. */
template<typename T> bool notStrPos(const T &a, std::string name) {
	if (a <= T(0)) {
		std::cerr << "Error: " << name << " should be strictly positive."
		          << std::endl;
		return true;
	}
	return false;
}

#endif // ACTIVEBROWNIAN_SIMUL_H_
