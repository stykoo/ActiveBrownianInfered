#ifndef ACTIVEBROWNIAN_BOXES_H_
#define ACTIVEBROWNIAN_BOXES_H_

#include <vector>
#include <array>

/* Class for the boxes in which particles are.  */
class Boxes {
	public:
		Boxes(const double Lx, const double Ly, const long n_parts,
			  const double size);
		//! Update according to the positions
		void update(const std::vector<double> &pos_x,
				    const std::vector<double> &pos_y);

		long getNBoxes() const {
			return n_boxes;
		}
		const auto& getNbrsPos() const {
			return nbrs_pos;
		}
		const auto& getPartsOfBox() const {
			return parts_of_box;
		}

	private:
		inline long ind(long i, long j) {
			return i * n_boxes_y + j;
		}

		const long n_boxes_x; //!< Number of boxes in the x direction
		const long n_boxes_y; //!< Number of boxes in the y direction
		const long n_boxes; //!< Total number of boxes
		const double len_box_x; //!< Length of a box in x direction 
		const double len_box_y; //!< Length of a box in y direction
		const long n_parts; //!< Number of particles
		//! Neighboring boxes of a given box along the 'positive' directions
		std::vector<std::array<long, 4>> nbrs_pos; 
		//!< Particles in a given box
		std::vector< std::vector<long> > parts_of_box;
};

#endif // ACTIVEBROWNIAN_BOXES_H_
