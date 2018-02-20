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

/*!
 * \brief Constructor for visualization
 *
 * \param state Pointer to the state of the system
 * \param len Length of the box
 * \param n_parts Number of particles
 */
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

    sf::CircleShape circle(scale);

    window.setFramerateLimit(FPS);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::White);

		for (long i = 0 ; i < n_parts ; ++i) {
			float x = state->getPos(i)[0] * scale; 
			pbc(x, (float) windowSize);
			float y = state->getPos(i)[1] * scale; 
			pbc(y, (float) windowSize);
			circle.setPosition(x, y);
			circle.setFillColor(colorFromAngle(state->getAngle(i)));
            window.draw(circle);
		}
        window.display();
    }
}

/*!
 * \brief Associate a color to an angle
 *
 * Use HSV representation with H = angle in degrees, S = 1, V = 1
 * 
 * \param angle Angle between 0 and 2pi
 * \return SFML color
 */
sf::Color colorFromAngle(const double angle) {
	const int hue = (int) (angle * 180. / M_PI);
	const float sat = 1.0;
	const float val = 1.0;

	int h = hue/60;
	float f = float(hue)/60-h;
	float p = val*(1.f-sat);
	float q = val*(1.f-sat*f);
	float t = val*(1.f-sat*(1-f));

	switch(h)
	{
		default:
		case 0:
		case 6:
			return sf::Color(val*255, t*255, p*255);
		case 1:
			return sf::Color(q*255, val*255, p*255);
		case 2:
			return sf::Color(p*255, val*255, t*255);
		case 3:
		 	return sf::Color(p*255, q*255, val*255);
		case 4:
			return sf::Color(t*255, p*255, val*255);
		case 5:
			return sf::Color(val*255, p*255, q*255);
	}
}
