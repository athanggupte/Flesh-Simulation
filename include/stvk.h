#pragma once

#include <Eigen/Dense>

namespace flesh {

	struct StVKMaterial
	{
		double lambda; // Lame's first parameter
		double mu; // Lame's second parameter

		void compute_piola_kirchhoff_stress(
			Eigen::Matrix3d const& F,
			Eigen::Matrix3d& P) const;

		auto get_piola_kirchhoff_fn() const
		{
			return [this](Eigen::Matrix3d const& F, Eigen::Matrix3d& P)
			{
				compute_piola_kirchhoff_stress(F, P);
			};
		}
	};

}
