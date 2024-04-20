#pragma once

namespace flesh
{
	/**
	 * \brief Computes the Lame coefficients from Young's modulus and Poisson ratio
	 * \param youngs_modulus Young's modulus
	 * \param poisson_ratio Poisson ratio
	 * \param lambda Lame's first parameter
	 * \param mu Lame's second parameter
	 */
	void lame_coefficients(
		double youngs_modulus,
		double poisson_ratio,
		double& lambda,
		double& mu)
	{
		mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
		lambda = youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
	}
}
