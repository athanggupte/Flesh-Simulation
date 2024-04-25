#pragma once

#include "base_material.h"
#include "piola_kirchhoff.h"

namespace flesh {

	struct StVKMaterial : public BaseMaterial
	{
		/**
		 * \brief Compute the first Piola-Kirchhoff stress tensor for St. Venant-Kirchhoff material
		 *
		 *     P(F) = F * ((2 * mu) * E + lambda * tr(E) * I)
		 *
		 *	   dP(F) = dF * ((2 * mu) * E + lambda * tr(E) * I) + F * ((2 * mu) * dE + lambda * tr(dE) * I)
		 */
		void compute_piola_kirchhoff_stress(
			Eigen::Matrix3d const& F,
			Eigen::Matrix3d& P) const override;

		void compute_piola_kirchhoff_stress_differential(
			Eigen::Matrix3d const& F,
			Eigen::Matrix3d const& dF,
			Eigen::Matrix3d& P,
			Eigen::Matrix9x9d& dPdF) const override;

		void precompute() override {}
	};

}
