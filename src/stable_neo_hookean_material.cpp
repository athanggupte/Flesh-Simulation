#include "stable_neo_hookean_material.h"

namespace flesh {

	void StableNeoHookeanMaterial::compute_piola_kirchhoff_stress(
		Eigen::Matrix3d const& F,
		Eigen::Matrix3d& P) const
	{
		// Compute the invariants of the deformation gradient
		double J = F.determinant();
		Eigen::Matrix3d C = F.transpose() * F;
		double IC = C.trace();

		Eigen::Matrix3d dJdF;
		dJdF.col(0) = F.col(1).cross(F.col(2));
		dJdF.col(1) = F.col(2).cross(F.col(0));
		dJdF.col(2) = F.col(0).cross(F.col(1));

		P = (mu * (1. - 1. / (IC + 1.))) * F + (lambda * (J - alpha)) * dJdF;
	}

	void StableNeoHookeanMaterial::compute_piola_kirchhoff_stress_differential(
		Eigen::Matrix3d const& F,
		Eigen::Matrix3d const& dF,
		Eigen::Matrix3d& P,
		Eigen::Matrix9x9d& dPdF) const
	{
	}

	void StableNeoHookeanMaterial::precompute()
	{
		alpha = 1. + mu / lambda - mu / lambda / 4.;
	}
}
