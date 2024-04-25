#pragma once
#include "base_material.h"
#include "piola_kirchhoff.h"

namespace flesh {

	struct StableNeoHookeanMaterial : public BaseMaterial
	{
		double alpha;

		/**
		 * \brief Compute the first Piola-Kirchhoff stress tensor for Stable Neo-Hookean material
		 */
		void compute_piola_kirchhoff_stress(
			Eigen::Matrix3d const& F,
			Eigen::Matrix3d& P) const override;

		void compute_piola_kirchhoff_stress_differential(
			Eigen::Matrix3d const& F,
			Eigen::Matrix3d const& dF,
			Eigen::Matrix3d& P,
			Eigen::Matrix9x9d& dPdF) const override;

		void precompute() override;
	};

}
