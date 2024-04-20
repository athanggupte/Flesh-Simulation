#include "..\include\shape_matrices.h"

namespace flesh {

	void shape_matrices(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXi const& T,
		Eigen::MatrixXd& D)
	{
		const Eigen::Index nF = T.rows(); // Number of tetrahedra

		D.resize(nF * 3, 3);

		for (Eigen::Index i = 0; i < nF; i++)
		{
			const Eigen::RowVector4i tet = T.row(i);

			Eigen::Vector3d v0 = X.row(tet(0));
			Eigen::Vector3d v1 = X.row(tet(1));
			Eigen::Vector3d v2 = X.row(tet(2));
			Eigen::Vector3d v3 = X.row(tet(3));

			Eigen::Matrix3d Di;
			Di.col(0) = v0 - v3;
			Di.col(1) = v1 - v3;
			Di.col(2) = v2 - v3;

			D.block<3, 3>(i * 3, 0) = Di;
		}
	}
}
