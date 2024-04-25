#pragma once

#include <Eigen/Dense>

namespace Eigen {
	using Matrix9x9d = Matrix<double, 9, 9>;
	using Matrix9x12d = Matrix<double, 9, 12>;
}

namespace flesh {
	/**
	 * \brief Function signature for computation of first Piola-Kirchhoff stress tensor
	 * \param F 3 x 3 deformation gradient
	 * \param P 3 x 3 first Piola-Kirchhoff stress tensor
	 */
	using piola_kirchhoff_strain_fn = std::function<void(
		Eigen::Matrix3d const& F,
		Eigen::Matrix3d& P)>;

	/**
	 * \brief Function signature for computation of first Piola-Kirchhoff stress tensor and its differential
	 * \param[in] F   3 x 3 deformation gradient
	 * \param[in] dF  3 x 3 gradient of the deformation gradient
	 * \param[out] P  3 x 3 first Piola-Kirchhoff stress tensor
	 * \param[out] dPdF 9 x 9 gradient of the first Piola-Kirchhoff stress tensor w.r.t. F
	 */
	using piola_kirchhoff_strain_differential_fn = std::function<void(
		Eigen::Matrix3d const& F,
		Eigen::Matrix3d const& dF,
		Eigen::Matrix3d& P,
		Eigen::Matrix9x9d& dPdF)>;

}
