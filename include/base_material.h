#pragma once
#include "piola_kirchhoff.h"

namespace flesh {

	struct BaseMaterial
	{
		double lambda; // Lame's first parameter
		double mu; // Lame's second parameter

		virtual ~BaseMaterial() = default;

		virtual void compute_piola_kirchhoff_stress(
			Eigen::Matrix3d const& F,
			Eigen::Matrix3d& P) const = 0;

		virtual void compute_piola_kirchhoff_stress_differential(
			Eigen::Matrix3d const& F,
			Eigen::Matrix3d const& dF,
			Eigen::Matrix3d& P,
			Eigen::Matrix9x9d& dPdF) const = 0;

		virtual void precompute() = 0;

		auto get_piola_kirchhoff_fn() const
		{
			return [this](Eigen::Matrix3d const& F, Eigen::Matrix3d& P)
			{
				compute_piola_kirchhoff_stress(F, P);
			};
		}

		auto get_piola_kirchhoff_differential_fn() const
		{
			return [this](Eigen::Matrix3d const& F, Eigen::Matrix3d const& dF, Eigen::Matrix3d& P, Eigen::Matrix9x9d& dP)
			{
				compute_piola_kirchhoff_stress_differential(F, dF, P, dP);
			};
		}
	};

}
