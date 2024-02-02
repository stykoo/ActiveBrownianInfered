#include <iostream>
#include <cmath>
#include "visu.h"

/* Constructor for visualization */
Visu::Visu(const State *state, const double Lx, const double Ly,
		   const long n_parts) :
	state(state), Lx(Lx), Ly(Ly), n_parts(n_parts) {
}

/*  Thread for visualization. */
void Visu::run() {
	sf::VideoMode mode = sf::VideoMode::getDesktopMode();
	const float windowSize = std::min(mode.width, mode.height) * 9 / 10;

    sf::RenderWindow window;
    window.create(sf::VideoMode(windowSize, (int) (windowSize * Ly / Lx)),
	              "Active Brownian Particles");

	float scale = windowSize / Ly;
	// We assume that the particles have diameter 1
    sf::CircleShape circle(scale / 2.0);
	// Line for showing the orientation of the particles
	sf::RectangleShape line(sf::Vector2f(scale / 2.0, 2.0));
	line.setFillColor(sf::Color::Black);

    window.setFramerateLimit(FPS);

	const std::vector<double> & pos_x = state->getPosX();
	const std::vector<double> & pos_y = state->getPosY();
	const std::vector<double> & angles = state->getAngles();

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::White);

		for (long i = 0 ; i < n_parts ; ++i) {
			float x = pos_x[i] * scale; 
			pbc(x, (float) windowSize);
			float y = pos_y[i] * scale; 
			pbc(y, (float) windowSize);

			// Angle is coded both as color and as arrow
			circle.setFillColor(colorFromAngle(angles[i]));
			line.setRotation(angles[i] * 180.0 / M_PI);

			// Draw multiple times if on the boundary
			int per_x = (x > windowSize - scale);
			int per_y = (y > windowSize - scale);

			for (int px = 0 ; px <= per_x ; ++px) {
				for (int py = 0 ; py <= per_y ; ++py) {
					circle.setPosition(x - px * windowSize,
					                   y - py * windowSize);
					window.draw(circle);
					line.setPosition(x - px * windowSize + scale / 2.0,
					                 y - py * windowSize + scale / 2.0);
					window.draw(line);
				}
			}
		}
        window.display();
    }
}

/* Associate a color to an angle
 * Use HSV representation with H = angle in degrees, S = 1, V = 1 */
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
