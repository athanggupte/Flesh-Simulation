#pragma once

#include <Eigen/Dense>

namespace flesh {
	/**
	 * \brief Computes a forward Euler step (using the semi-implicit Euler integration scheme)
	 * \param[in] X #X x 3 vertex positions matrix
	 * \param[in] V #X x 3 vertex velocities matrix
	 * \param[in] f #X x 3 forces matrix
	 * \param[in] h Time step
	 * \param[out] X_next #X x 3 next vertex positions matrix
	 * \param[out] V_next #X x 3 next vertex velocities matrix
	 */
	void forward_euler_step(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXd const& V,
		Eigen::MatrixXd const& f,
		double h,
		Eigen::MatrixXd& X_next,
		Eigen::MatrixXd& V_next);
}
