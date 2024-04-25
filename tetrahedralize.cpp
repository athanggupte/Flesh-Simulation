#include <igl/copyleft/tetgen/tetrahedralize.h>

void tetrahedralize(
	Eigen::MatrixXd const& V,
	Eigen::MatrixXi const& F,
	double max_tet_volume,
	Eigen::MatrixXd& TV,
	Eigen::MatrixXi& TT,
	Eigen::MatrixXi& TF)
{
	std::string switches = "pq1.25,21a" + std::to_string(max_tet_volume);

	igl::copyleft::tetgen::tetrahedralize(V,F, switches, TV, TT, TF);
}
