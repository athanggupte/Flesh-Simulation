#include "flesh.h"
#include "..\include\shape_matrices.h"

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
		piola_kirchhoff_fn const& P_fn,
		Eigen::MatrixXd& f)
	{
		const Eigen::Index nV = X.rows(); // Number of vertices
		const Eigen::Index nF = T.rows(); // Number of tetrahedra

		f.setZero(nV, 3); // Elastic forces

		Eigen::MatrixXd D; // #T * 3 x 3 per-tet deformed shape matrices
		shape_matrices(X, T, D);

		for (Eigen::Index i = 0; i < nF; i++)
		{
			const Eigen::RowVector4i tet = T.row(i);

			const Eigen::Matrix3d Bi = B.block<3, 3>(i * 3, 0);
			const Eigen::Matrix3d Di = D.block<3, 3>(i * 3, 0);
			
			const Eigen::Matrix3d Fi = Di * Bi; // Deformation gradient

			Eigen::Matrix3d P; // First Piola-Kirchhoff stress tensor
			P_fn(Fi, P);

			const Eigen::Matrix3d H = -W(i) * P * Bi.transpose(); // Elastic forces

			f.row(tet(0)) += H.col(0);
			f.row(tet(1)) += H.col(1);
			f.row(tet(2)) += H.col(2);
			f.row(tet(3)) -= H.col(0) + H.col(1) + H.col(2);
		}
	}
}
