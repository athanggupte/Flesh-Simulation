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
	 * \param[in]  piola_kirchhoff_strain Function to compute the first Piola-Kirchhoff stress tensor
	 * \param[out] f #V x 3 per-vertex elastic forces
	 */
	void flesh_elastic_forces(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXi const& T,
		Eigen::MatrixXd const& B,
		Eigen::VectorXd const& W,
		piola_kirchhoff_strain_fn const& piola_kirchhoff_strain,
		Eigen::MatrixXd& f);

	/**
	 * \brief Computes the elastic forces and their differentials for the flesh simulation
	 * \param[in]  X #X x 3 vertex positions matrix
	 * \param[in]  T #T x 4 tetrahedral face matrix
	 * \param[in]  V #X x 3 vertex velocities matrix
	 * \param[in]  B #T * 3 x 3 per-tet inverse reference shape matrices
	 * \param[in]  W #F x 1 per-tet volume
	 * \param[in]  piola_kirchhoff_strain_differential Function to compute the first Piola-Kirchhoff stress tensor and its differential
	 * \param[out] f #V x 3 per-vertex elastic forces
	 * \param[out] df #V x 3 x 3 per-vertex elastic force differentials
	 */
	void flesh_force_differentials(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXi const& T,
		Eigen::MatrixXd const& V,
		Eigen::MatrixXd const& B,
		Eigen::VectorXd const& W,
		piola_kirchhoff_strain_differential_fn const& piola_kirchhoff_strain_differential,
		Eigen::MatrixXd& f,
		Eigen::MatrixXd& df);

}
