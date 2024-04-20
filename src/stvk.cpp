#include "stvk.h"

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
}
