#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace flesh {
	/**
	 * \brief Generates the shape matrix per tetrahedron of a mesh
	 * \param[in]  X #X x 4 vertex position matrix
	 * \param[in]  T #T x 4 tetrahedral face matrix
	 * \param[out] D #T * 3 x 3 per-tet shape matrices
	 */
	void shape_matrices(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXi const& T,
		Eigen::MatrixXd & D);

}
