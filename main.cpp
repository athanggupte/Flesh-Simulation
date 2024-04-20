#include <igl/opengl/glfw/Viewer.h>

#include <igl/copyleft/tetgen/tetrahedralize.h>

#include "flesh.h"
#include "forward_euler_step.h"
#include "lame_coefficients.h"
#include "stvk.h"


int main(int argc, char *argv[])
{
	// Inline mesh of a cube
	const Eigen::MatrixXd V = (Eigen::MatrixXd(8,3)<<
		0.0,0.0,0.0,
		0.0,0.0,1.0,
		0.0,1.0,0.0,
		0.0,1.0,1.0,
		1.0,0.0,0.0,
		1.0,0.0,1.0,
		1.0,1.0,0.0,
		1.0,1.0,1.0).finished();
	const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
		0,6,4,
		0,2,6,
		0,3,2,
		0,1,3,
		2,7,6,
		2,3,7,
		4,6,7,
		4,7,5,
		0,4,5,
		0,5,1,
		1,5,7,
		1,7,3).finished();

	Eigen::MatrixXd TV;
	Eigen::MatrixXi TT, TF;

	struct {
		Eigen::MatrixXd X; // Vertex positions
		Eigen::MatrixXi T; // Tetrahedral faces
		Eigen::MatrixXd V; // Vertex velocities
		Eigen::MatrixXd B; // Per-tet inverse reference shape matrices
		Eigen::VectorXd W; // Per-tet volume
	} simdata;
	flesh::StVKMaterial stvk;
	flesh::lame_coefficients(1e0, 0.49, stvk.lambda, stvk.mu);

	igl::copyleft::tetgen::tetrahedralize(V,F, "pq1.414Y", TV, TT, TF);

	simdata.X = TV;
	simdata.T = TT;
	simdata.V.setZero(simdata.X.rows(), simdata.X.cols());
	flesh::flesh_precompute(simdata.X, simdata.T, simdata.B, simdata.W);

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(simdata.X, TF);
	viewer.data().set_face_based(true);
	viewer.data().invert_normals = true;

	const auto deformed = [&simdata]()
	{
		simdata.X.row(1) -= Eigen::RowVector3d(0.25, 0.25, 0.25);
		return true;
	}();

	const auto simulation_update = [&simdata,&stvk](double dt)
	{
		Eigen::MatrixXd f_elastic;
		flesh::flesh_elastic_forces(simdata.X, simdata.T, simdata.B, simdata.W, stvk.get_piola_kirchhoff_fn(), f_elastic);

		Eigen::MatrixXd X, V;
		X.resizeLike(simdata.X);
		V.resizeLike(simdata.V);

		static double timer = 0.0;
		constexpr double step = 0.0001;
		timer += std::min(dt, 0.01);

		while (timer > step + FLT_EPSILON)
		{
			flesh::forward_euler_step(simdata.X, simdata.V, f_elastic, step, X, V);
			timer -= step;
			simdata.X = X;
			simdata.V = V;
		}
	};

	viewer.callback_pre_draw = [&simulation_update,&TV,&simdata](igl::opengl::glfw::Viewer &viewer)->bool
	{
		static double t_start = igl::get_seconds();
		double t = igl::get_seconds();
		double dt = t - t_start;
		t_start = t;

		// Animate the deformation
		simulation_update(dt);

		TV = simdata.X;

		viewer.data().set_vertices(TV);
		viewer.data().compute_normals();

		Eigen::MatrixXd P = TV.topRows(1);
		viewer.data().set_points(P, Eigen::RowVector3d(1.00, 0.34, 0.12));

		return false;
	};

	viewer.core().is_animating = true;

	viewer.launch();
}
