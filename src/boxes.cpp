#include <iostream>
#include <cmath>
#include "boxes.h"

/* Constructor of boxes. */
Boxes::Boxes(const double Lx, const double Ly, const long _n_parts,
	         const double size) :
		n_boxes_x(std::floor(Lx / size)),
		n_boxes_y(std::floor(Ly / size)),
		n_boxes(n_boxes_x * n_boxes_y),
		len_box_x(Lx / n_boxes_x),
		len_box_y(Ly / n_boxes_y),
		n_parts(_n_parts) {
	nbrs_pos.resize(n_boxes);

	for (long i = 0 ; i < n_boxes_x ; ++i) {
		for (long j = 0 ; j < n_boxes_y ; ++j) {
			long ip = (i + 1) % n_boxes_x;
			long jp = (j + 1) % n_boxes_y;
			long im = (i - 1 + n_boxes_x) % n_boxes_x;
			nbrs_pos[ind(i, j)][0] = ind(ip, j);
			nbrs_pos[ind(i, j)][1] = ind(i, jp);
			nbrs_pos[ind(i, j)][2] = ind(ip, jp);
			nbrs_pos[ind(i, j)][3] = ind(im, jp);
		}
	}

    /*for (long k = 0 ; k < n_boxes ; ++k) {
		std::cout << "[" << k << "] ";
		for (long i : nbrs_pos[k]) {
			std::cout << i << " ";
		}
		std::cout << std::endl;
	}*/

	parts_of_box.resize(n_boxes);
}

/* Classify the particles at given positions in the boxes.
 * Warning: the positions should be consistent with the Boxes object!  */
void Boxes::update(const std::vector<double> &pos_x,
		           const std::vector<double> &pos_y) {
    for (long i = 0 ; i < n_boxes ; ++i) {
		parts_of_box[i].clear();
	}

    for (long i = 0 ; i < n_parts ; ++i) {
		// Round towards 0
		long bx = (long) (pos_x[i] / len_box_x);
		long by = (long) (pos_y[i] / len_box_y);
		parts_of_box[ind(bx, by)].push_back(i);
	}
}
