#pragma once

#include <Eigen/Dense>

namespace flesh {
	/**
	 * \brief Function signature for the first Piola-Kirchhoff stress tensor computation
	 * \param F 3 x 3 deformation gradient
	 * \param P 3 x 3 first Piola-Kirchhoff stress tensor
	 */
	using piola_kirchhoff_fn = std::function<void(
		Eigen::Matrix3d const& F,
		Eigen::Matrix3d& P)>;

}
