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
 * \brief Visualization of the system
 *
 * The system is visualized using the SFML library.
*/

#include "visu.h"
#include <iostream>

Visu::Visu(const State *state, const double len, const long n_parts) :
	state(state), len(len), n_parts(n_parts), scale(windowSize / len) {
}

/*!
 * \brief Thread for visualization.
 *
 * Open a window, draw the particles and update their
 * positions at a certain number of FPS while the simulation is runing.
 */
void Visu::run() {
    sf::RenderWindow window;
    window.create(sf::VideoMode(windowSize, windowSize),
	              "Active Brownian Particles");

    sf::CircleShape circle(Visu::circleRad);
	circle.setFillColor(sf::Color::Green);

    window.setFramerateLimit(Visu::FPS);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::Black);

		for (long i = 0 ; i < n_parts ; ++i) {
			float x = state->getPos(i)[0] * scale; 
			pbc(x, (float) windowSize);
			float y = state->getPos(i)[1] * scale; 
			pbc(y, (float) windowSize);
			circle.setPosition(x, y);
            window.draw(circle);
		}
        window.display();
    }
}
