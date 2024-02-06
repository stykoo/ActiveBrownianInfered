#ifndef ACTIVEBROWNIAN_STATE_H_
#define ACTIVEBROWNIAN_STATE_H_

#include <vector>
#include <array>
#include "infered.h"
#include "boxes.h"

#ifdef USE_MKL
	#include "mkl.h"
	#include "mkl_vsl.h"
#else
	#include <random>
#endif

// 2^(1/6)
#define TWOONESIXTH 1.12246204830937298143 

/*!  Class for the state of the system */
class State {
	public:
		//! Constructor of State
		State(const double _Lx, const double _Ly, const long _n_parts,
			  const double _diam, const double _trans_dif,
			  const double _rot_dif, const double _activity, const double _dt,
			  const double _pot_strength, Infered &_infered);
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
			return positions[2];
		}

		void store(std::vector<double> &out) const; //!< Store positions
		void dump() const; //!< Dump the positions and orientations

	private:
		void calcInternalForces(); //!< Compute internal forces
		void calcSoftForceIJ(const long i, const long j);
		void calcInferedForceIJ(const long i, const long j);
		void enforcePBC(); //!< Enforce periodic boundary conditions

		const std::array<double, 3> lengths;
		const long n_parts; //!< Number of particles
		double diam; //!< Particle diameter
		const double activity; //!< Activity
		const double dt; //!< Timestep
		const double pot_strength; //!< Strength of the repulsive potential
		Infered &infered; //!< Structure for inference
		Boxes boxes; //!< Structure for boxes

#ifdef USE_MKL
		std::array<double, 3> stddevs;
		VSLStreamStatePtr stream;
		std::array<std::vector<double>, 3> aux;
#else
		std::mt19937 rng; //!< Random number generator
		//! Gaussian noise for translational diffusivity
		std::normal_distribution<double> noiseTemp;
		//! Gaussian noise for angle
		std::normal_distribution<double> noiseAngle;
#endif

		//! Positions of the particles and angles
		std::array<std::vector<double>, 3> positions;
  		//!< Internal forces
		std::array<std::vector<double>, 3> forces;
};

#ifdef USE_MKL
void pbcMKL(std::vector<double> &v, const double L, std::vector<double> &aux,
	        const long N);
void pbcSymMKL(std::vector<double> &v, const double L,
		       std::vector<double> &aux, const long N);
#endif

/* Periodic boundary conditions on a segment */
template<typename T>
void pbc(T &x, const T L){
	x -= L * std::floor(x / L);
}

/* Periodic boundary conditions on a segment (symmetric) */
template<typename T>
inline void pbcSym(T &x, const T L) {
	x -= L * std::round(x / L);
}

/*
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
*/

#endif // ACTIVEBROWNIAN_STATE_H_
