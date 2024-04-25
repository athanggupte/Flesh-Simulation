#include "flesh.h"

#include <iostream>

#include "shape_matrices.h"

#include <omp.h>

namespace flesh {

	void flesh_precompute(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXi const& T,
		Eigen::MatrixXd& B,
		Eigen::VectorXd& W)
	{
		const Eigen::Index nF = T.rows(); // Number of tetrahedra

		Eigen::MatrixXd D;
		shape_matrices(X, T, D);

		B.resize(nF * 3, 3);
		W.resize(nF);

		for (Eigen::Index i = 0; i < nF; i++)
		{
			const Eigen::Matrix3d Di = D.block<3, 3>(i * 3, 0);
			W(i) = abs(Di.determinant()) / 6.0; // Volume of the tetrahedron
			B.block<3, 3>(i * 3, 0) = Di.inverse(); // Inverse of the reference shape matrix
		}
	}

	void flesh_elastic_forces(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXi const& T,
		Eigen::MatrixXd const& B,
		Eigen::VectorXd const& W,
		piola_kirchhoff_strain_fn const& piola_kirchhoff_strain,
		Eigen::MatrixXd& f)
	{
		const Eigen::Index nV = X.rows(); // Number of vertices
		const Eigen::Index nF = T.rows(); // Number of tetrahedra

		f.setZero(nV, 3); // Elastic forces

		Eigen::MatrixXd D; // #T * 3 x 3 per-tet deformed shape matrices
		shape_matrices(X, T, D);

		Eigen::MatrixXd H(nF * 3, 3); // Elastic forces per tetrahedron (3 vertices per tetrahedron)

		#pragma omp parallel for
		for (Eigen::Index i = 0; i < nF; i++)
		{
			const Eigen::Matrix3d Bi = B.block<3, 3>(i * 3, 0);
			const Eigen::Matrix3d Di = D.block<3, 3>(i * 3, 0);
			
			const Eigen::Matrix3d Fi = Di * Bi; // Deformation gradient

			Eigen::Matrix3d P; // First Piola-Kirchhoff stress tensor
			piola_kirchhoff_strain(Fi, P);

			H.block<3, 3>(i * 3, 0) = -W(i) * P * Bi.transpose(); // Elastic forces
		}

		for (Eigen::Index i = 0; i < nF; i++)
		{
			const Eigen::RowVector4i tet = T.row(i);
			const Eigen::Matrix3d Hi = H.block<3, 3>(i * 3, 0);
			f.row(tet(0)) += Hi.col(0);
			f.row(tet(1)) += Hi.col(1);
			f.row(tet(2)) += Hi.col(2);
			f.row(tet(3)) -= Hi.col(0) + Hi.col(1) + Hi.col(2);
		}
	}

	void flesh_force_differentials(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXi const& T,
		Eigen::MatrixXd const& V,
		Eigen::MatrixXd const& B,
		Eigen::VectorXd const& W,
		piola_kirchhoff_strain_differential_fn const& piola_kirchhoff_strain_differential,
		Eigen::MatrixXd& f,
		Eigen::MatrixXd& df)
	{
		
	}
}
