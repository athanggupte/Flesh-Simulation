#pragma once

#include <Eigen/Dense>

#include "piola_kirchhoff.h"

namespace flesh {
	/**
	 * \brief Precomputes initial matrices for the flesh simulation
	 * \param[in]  X #X x 3 vertex positions matrix
	 * \param[in]  T #T x 4 tetrahedral face matrix
	 * \param[out] B #T * 3 x 3 per-tet inverse reference shape matrices
	 * \param[out] W #F x 1 per-tet volume
	 */
	void flesh_precompute(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXi const& T,
		Eigen::MatrixXd& B,
		Eigen::VectorXd& W);

	/**
	 * \brief Computes the elastic forces for the flesh simulation
	 * \param[in]  X #X x 3 vertex positions matrix
	 * \param[in]  T #T x 4 tetrahedral face matrix
	 * \param[in]  B #T * 3 x 3 per-tet inverse reference shape matrices
	 * \param[in]  W #F x 1 per-tet volume
	 * \param[in]  P_fn Function to compute the first Piola-Kirchhoff stress tensor
	 * \param[out] f #V x 3 per-vertex elastic forces
	 */
	void flesh_elastic_forces(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXi const& T,
		Eigen::MatrixXd const& B,
		Eigen::VectorXd const& W,
		piola_kirchhoff_fn const& P_fn,
		Eigen::MatrixXd& f);

}
