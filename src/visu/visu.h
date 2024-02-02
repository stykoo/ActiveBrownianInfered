#ifndef ACTIVEBROWNIAN_VISU_H_
#define ACTIVEBROWNIAN_VISU_H_

#include <SFML/Graphics.hpp>
#include "../state.h"

class Visu {
	public:
		Visu(const State *state, const double Lx, const double Ly,
			 const long n_parts);
		void run();

	private:
		const int FPS = 24; //!< Number of frames per second

		const State *state; //!< Pointer to the state of the system
		const double Lx, Ly; //!< Lengths of the box
		const long n_parts; //!< Number of particles
};

sf::Color colorFromAngle(const double angle);

#endif // ACTIVEBROWNIAN_VISU_H_
