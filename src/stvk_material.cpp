#include "..\include\stvk_material.h"

namespace flesh {

	void StVKMaterial::compute_piola_kirchhoff_stress(
		Eigen::Matrix3d const& F,
		Eigen::Matrix3d& P) const
	{
		Eigen::Matrix3d F_T = F.transpose();

		Eigen::Matrix3d E = 0.5 * (F_T * F - Eigen::Matrix3d::Identity()); // Green strain tensor

		Eigen::Matrix3d lambda_trE = lambda * E.trace() * Eigen::Matrix3d::Identity();

		P = F * ((2.0 * mu) * E + lambda_trE);
	}

	void StVKMaterial::compute_piola_kirchhoff_stress_differential(
		Eigen::Matrix3d const& F,
		Eigen::Matrix3d const& dF,
		Eigen::Matrix3d& P,
		Eigen::Matrix9x9d& dPdF) const
	{
		Eigen::Matrix3d F_T = F.transpose();
		Eigen::Matrix3d dF_T = dF.transpose();

		Eigen::Matrix3d E = 0.5 * (F_T * F - Eigen::Matrix3d::Identity()); // Green strain tensor
		Eigen::Matrix3d dE = 0.5 * (dF_T * F + F_T * dF);

		Eigen::Matrix3d lambda_trE = lambda * E.trace() * Eigen::Matrix3d::Identity();
		Eigen::Matrix3d lambda_trdE = lambda * dE.trace() * Eigen::Matrix3d::Identity();

		P = F * ((2.0 * mu) * E + lambda_trE);
		//dPdF = dF * ((2.0 * mu) * E + lambda_trE) + F * ((2.0 * mu) * dE + lambda_trdE);// TODO: This is not correct. This is dP and not dPdF
	}
}
