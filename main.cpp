#include <igl/opengl/glfw/Viewer.h>

#include "flesh.h"
#include "forward_euler_step.h"
#include "lame_coefficients.h"
#include "include\stvk_material.h"
#include "getopt.h"
#include "stable_neo_hookean_material.h"

void tetrahedralize(
	Eigen::MatrixXd const& V,
	Eigen::MatrixXi const& F,
	double max_tet_volume,
	Eigen::MatrixXd& TV,
	Eigen::MatrixXi& TT,
	Eigen::MatrixXi& TF);

void make_cube_bar(int side, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
	V.resize((side + 1) * 4, 3);
	F.resize(side * 8 + 4, 3);

	// makes a bar side x 1 x 1 cubes
	for (int s = 0; s < side + 1; s++)
	{
		V.row(s * 4 + 0) << s, 0, 0;
		V.row(s * 4 + 1) << s, 0, 1;
		V.row(s * 4 + 2) << s, 1, 1;
		V.row(s * 4 + 3) << s, 1, 0;
	}

	F.row(0) << 0, 1, 2;
	F.row(1) << 0, 2, 3;

	for (int c = 0; c < side; c++)
	{
		int v[8];
		for (int i = 0; i < 8; i++)
		{
			v[i] = c * 4 + i;
		}
		
		F.row(c * 8 + 2) << v[1], v[5], v[6];
		F.row(c * 8 + 3) << v[1], v[6], v[2];
		F.row(c * 8 + 4) << v[0], v[5], v[4];
		F.row(c * 8 + 5) << v[0], v[5], v[1];
		F.row(c * 8 + 6) << v[0], v[3], v[4];
		F.row(c * 8 + 7) << v[3], v[7], v[4];
		F.row(c * 8 + 8) << v[3], v[6], v[7];
		F.row(c * 8 + 9) << v[3], v[2], v[6];
	}

	F.row(side * 8 + 2) << side * 4 + 0, side * 4 + 1, side * 4 + 2;
	F.row(side * 8 + 3) << side * 4 + 0, side * 4 + 2, side * 4 + 3;
}

enum FixtureType
{
	eFixtureType_None,
	eFixtureType_OneEnd,
	eFixtureType_BothEnds
};

enum MaterialType
{
	eMaterialType_StVK,
	eMaterialType_StableNeoHookean
};

void help()
{
	std::cout << R"(Usage:    driver [options...]

    Options:
  -n <frames>         Number of frames to write to file
  -w                  Write simulation to file
  -p <file>           Playback simulation from file
  -f <file>           Model file to load (default: cube_bar)
  -Y <youngs modulus> Young's modulus (default: 1e6)
  -P <poisson ratio>  Poisson ratio (default: 0.49)
  -m <material type>  Material type [stvk, stneo] (default: stneo)
)";
}

int main(int argc, char *argv[])
{
	static bool write_to_file = false;
	static bool playback_from_file = false;
	static int  max_frames = -1;
	static bool is_simulating = false;
	static bool is_looping = true;
	static std::string playback_file = "sim.dat";
	static std::string model_file;
	static FixtureType fixture_type = eFixtureType_BothEnds;
	static MaterialType material_type = eMaterialType_StableNeoHookean;
	static double fixture_thickness = 0.2;
	static double youngs_modulus = 1e6;
	static double poisson_ratio = 0.49;
	static bool gravity = true;

	static int frameNumber = 0;

	if (argc > 1)
	{
		int c;
		while ((c = getopt(argc, argv, "n:wp:f:Y:P:m:")) != -1)
		{
			switch (c)
			{
			case 'n':
				max_frames = atoi(optarg);
				break;
			case 'w':
				write_to_file = true;
				break;
			case 'p':
				playback_from_file = true;
				playback_file = optarg;
				break;
			case 'f':
				model_file = optarg;
				break;
			case 'Y':
				youngs_modulus = atof(optarg);
				break;
			case 'P':
				poisson_ratio = atof(optarg);
				break;
			case 'm':
				if (strcmp(optarg, "stvk") == 0)
				{
					material_type = eMaterialType_StVK;
				}
				else if (strcmp(optarg, "stneo") == 0)
				{
					material_type = eMaterialType_StableNeoHookean;
				}
				break;
			case '?':
				help();
				exit(1);
			}
		}
	}

	if (write_to_file && playback_from_file)
	{
		std::cerr << "Cannot write to file and playback from file at the same time" << std::endl;
		exit(1);
	}

	if (write_to_file && max_frames < 0)
	{
		std::cerr << "Must specify number of frames to write with -n" << std::endl;
		exit(1);
	}

	Eigen::MatrixXd TV;
	Eigen::MatrixXi TT, TF;

	if (playback_from_file)
	{
		igl::deserialize(max_frames, "max_frames", playback_file);
		igl::deserialize(TV, "TV", playback_file);
		igl::deserialize(TT, "TT", playback_file);
		igl::deserialize(TF, "TF", playback_file);
		
		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(TV, TF);
		viewer.data().set_face_based(true);
		viewer.data().invert_normals = true;
		viewer.core().is_animating = true;

		viewer.callback_pre_draw = [&TV](igl::opengl::glfw::Viewer &viewer)->bool
		{
			char buf[7] = {0};
			sprintf(buf, "TV%04d", frameNumber++);

			igl::deserialize(TV, std::string(buf), playback_file);

			viewer.data().set_vertices(TV);
			viewer.data().compute_normals();
			return false;
		};

		viewer.callback_post_draw = [&TV](igl::opengl::glfw::Viewer &viewer)->bool
		{
			if (frameNumber >= max_frames)
			{
				if (is_looping)
				{
					frameNumber = 0;
				}
				else
				{
					exit(0);
				}
			}

			return false;
		};

		viewer.launch();

		return 0;
	}
	
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	if (model_file.empty())
	{
		make_cube_bar(8, V, F);
		Eigen::Vector3d side_len(4., 4., 4.);
		V.col(0) *= side_len(0);
		V.col(1) *= side_len(1);
		V.col(2) *= side_len(2);
	}
	else
	{
		if (!igl::readOBJ(model_file, V, F))
		{
			std::cerr << "Failed to read model file" << std::endl;
			exit(1);
		}
	}

	tetrahedralize(V, F, 0.5, TV, TT, TF);

	struct {
		Eigen::MatrixXd X; // Vertex positions
		Eigen::MatrixXi T; // Tetrahedral faces
		Eigen::MatrixXd V; // Vertex velocities
		Eigen::MatrixXd B; // Per-tet inverse reference shape matrices
		Eigen::VectorXd W; // Per-tet volume
	} simdata;

	flesh::BaseMaterial * material = nullptr;
	switch (material_type)
	{
	case eMaterialType_StVK:
		std::cout << "Using StVK material" << '\n';
		material = new flesh::StVKMaterial();
		break;
	case eMaterialType_StableNeoHookean:
		std::cout << "Using Stable Neo-Hookean material" << '\n';
		material = new flesh::StableNeoHookeanMaterial();
		break;
	}

	flesh::lame_coefficients(youngs_modulus, poisson_ratio, material->lambda, material->mu);
	material->precompute();

	simdata.X = TV;
	simdata.T = TT;
	simdata.V.setZero(simdata.X.rows(), simdata.X.cols());
	flesh::flesh_precompute(simdata.X, simdata.T, simdata.B, simdata.W);

	if (write_to_file)
	{
		igl::serialize(max_frames, "max_frames", "sim.dat");
		igl::serialize(TV, "TV", "sim.dat");
		igl::serialize(TT, "TT", "sim.dat");
		igl::serialize(TF, "TF", "sim.dat");
	}

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(TV, TF);
	viewer.data().set_face_based(true);
	viewer.data().invert_normals = true;
	viewer.core().is_animating = true;

	const auto deformed = [&simdata]()
	{
		//simdata.X.row(simdata.X.rows()-1) += Eigen::RowVector3d(0, 0.1, 0);
		return true;
	}();

	Eigen::MatrixXd f_external;
	f_external.resizeLike(simdata.X);
	for (int i = 0; i < f_external.rows(); i++)
	{
		f_external.row(i) = Eigen::Vector3d(0, -9.8, 0);
	}

	static std::vector<Eigen::Index> edge_face_vertices;
	double max_x = V.colwise().maxCoeff()(0);
	double min_x = V.colwise().minCoeff()(0);
	static int left_fixed_count = 0, right_fixed_count = 0;

	// Fix left end
	if (fixture_type == eFixtureType_OneEnd || fixture_type == eFixtureType_BothEnds)
	{
		for (int i = 0; i < simdata.X.rows(); i++)
		{
			if (simdata.X(i, 0) <= min_x + fixture_thickness + DBL_EPSILON)
			{
				edge_face_vertices.push_back(i);
			}
		}
	}
	left_fixed_count = edge_face_vertices.size();

	// Fix right end
	if (fixture_type == eFixtureType_BothEnds)
	{
		for (int i = 0; i < simdata.X.rows(); i++)
		{
			if (simdata.X(i, 0) >= max_x - fixture_thickness - DBL_EPSILON)
			{
				edge_face_vertices.push_back(i);
			}
		}
	}
	right_fixed_count = edge_face_vertices.size() - left_fixed_count;

	static Eigen::Matrix<double,Eigen::Dynamic,3> edge_verts_start(edge_face_vertices.size(), 3);
	for (int i = 0; i < edge_face_vertices.size(); i++)
	{
		edge_verts_start.row(i) = simdata.X.row(edge_face_vertices[i]);
	}
	static Eigen::RowVector3d edge_center = edge_verts_start.colwise().mean();
	static Eigen::Matrix<double,Eigen::Dynamic,3> edge_radial = edge_verts_start.rowwise() - edge_center;

	static double angle = 0;

	const auto simulation_update = [&simdata,&material,&TV,&f_external](double dt)
	{
		const double deltaTime = std::min(dt, 1./60.);

		static double time = 0;
		time += deltaTime;

		static double timer = 0.0;
		constexpr double step = 0.0001;
		timer += deltaTime;

		angle = std::min(time * igl::PI / 16, igl::PI);

		while (timer > step + FLT_EPSILON)
		{
		#if TWIST_TEST
			Eigen::Matrix<double,Eigen::Dynamic,3> edge_verts(edge_face_vertices.size(), 3);
			for (int i = 0; i < edge_face_vertices.size(); i++)
			{
				edge_verts.row(i) = simdata.X.row(edge_face_vertices[i]);
			}
			Eigen::RowVector3d edge_center = edge_verts.colwise().mean();
			Eigen::RowVector3d edge_normal = edge_radial.row(0).cross(edge_radial.row(1)).normalized();
			Eigen::Matrix<double,Eigen::Dynamic,3> edge_tangent = edge_radial.rowwise().cross(edge_normal).normalized();
			for (int i = 0; i < edge_face_vertices.size(); i++)
			{
				//auto rotated_radial_edge = edge_radial.row(i) + edge_tangent.row(i) * step;
				// edge_verts.row(i) = rotated_radial_edge.normalized() * edge_radial.row(i).norm() + edge_center;
				Eigen::RowVector3d rotated_radial_edge = Eigen::AngleAxisd(angle, edge_normal).toRotationMatrix() * edge_radial.row(i).transpose();
				edge_verts.row(i) = rotated_radial_edge + edge_center;

				TV.row(edge_face_vertices[i]) = edge_verts.row(i);
			}
		#endif

		#if SHAKE_TEST || 1
			for (int i = 0; i < edge_face_vertices.size(); i++)
			{
				TV.row(edge_face_vertices[i]) = edge_verts_start.row(i) + Eigen::RowVector3d(0, 0, 2 * sin(time * 2));
			}

		#endif

			for (int i = 0; i < edge_face_vertices.size(); i++)
			{
				simdata.X.row(edge_face_vertices[i]) = TV.row(edge_face_vertices[i]);
			}

			if (gravity)
			{
				for (int i = 0; i < f_external.rows(); i++)
				{
					f_external.row(i) = Eigen::Vector3d(0, -9.8, 0);
				}
			}
			else
			{
				f_external.setZero();
			}

			Eigen::MatrixXd f_elastic;
			//f_elastic.setZero(simdata.X.rows(), 3);
			flesh::flesh_elastic_forces(simdata.X, simdata.T, simdata.B, simdata.W, material->get_piola_kirchhoff_fn(), f_elastic);

			// Eigen::VectorXd f_elastic_norm = f_elastic.rowwise().norm();
			double f_elastic_norm_max = f_elastic.rowwise().norm().maxCoeff<Eigen::PropagateNaN>();

			if (f_elastic.hasNaN())
			{
			#if _DEBUG
				__debugbreak();
			#endif
				std::cerr << "Forces are NaN!" << '\n';
				exit(1);
			}
			if (f_elastic_norm_max > 1e300)
			{
			#if _DEBUG
				__debugbreak();
			#endif
				std::cerr << "Forces too large!" << '\n';
				exit(1);
			}


			Eigen::MatrixXd X, V;
			X.resizeLike(simdata.X);
			V.resizeLike(simdata.V);

			/*for (int i = 0; i < edge_face_vertices.size(); i++)
			{
				f_external.row(edge_face_vertices[i]) = edge_tangent.row(i) * 10000.;
			}*/

			Eigen::MatrixXd f_damping = -simdata.V * 0.05;

			flesh::forward_euler_step(simdata.X, simdata.V, f_elastic + f_external + f_damping, step, X, V);
			timer -= step;
			simdata.X = X;
			simdata.V = V;

			for (int i = 0; i < edge_face_vertices.size(); i++)
			{
				simdata.X.row(edge_face_vertices[i]) = TV.row(edge_face_vertices[i]);
				simdata.V.row(edge_face_vertices[i]) = Eigen::RowVector3d(0, 0, 0);
			}
		}
	};

	viewer.callback_pre_draw = [&simulation_update,&TV,&simdata](igl::opengl::glfw::Viewer &viewer)->bool
	{
		static double t_start = igl::get_seconds();
		double t = igl::get_seconds();
		double dt = t - t_start;
		t_start = t;

		// Animate the deformation
		if (is_simulating)
		{
			simulation_update(dt);
			TV = simdata.X;
		}

		viewer.data().set_vertices(TV);
		viewer.data().compute_normals();

		Eigen::Matrix<double,Eigen::Dynamic,3> edge_verts(edge_face_vertices.size(), 3);
		for (int i = 0; i < edge_face_vertices.size(); i++)
		{
			edge_verts.row(i) = simdata.X.row(edge_face_vertices[i]);
		}
		viewer.data().set_points(edge_verts, Eigen::RowVector3d(0.87, 0.35, 0.21));

		return false;
	};

	viewer.callback_post_draw = [&TV](igl::opengl::glfw::Viewer &viewer)->bool
	{
		if (is_simulating && write_to_file)
		{
			char buf[7] = {0};
			sprintf(buf, "TV%04d", frameNumber++);
			igl::serialize(TV, std::string(buf), "sim.dat");

			if (frameNumber % 10 == 0)
			{
				printf("Frames written: %d\n", frameNumber);
			}
		}

		if (max_frames > 0 && frameNumber >= max_frames)
		{
			exit(0);
		}

		return false;
	};

	viewer.callback_key_down = [](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod)->bool
	{
		switch (key)
		{
		case GLFW_KEY_SPACE:
			is_simulating = !is_simulating;
			break;
		case GLFW_KEY_R:
			is_looping = !is_looping;
			break;
		case GLFW_KEY_G:
			gravity = !gravity;
			std::cout << "Gravity: " << (gravity ? "On" : "Off") << '\n';
			break;
		}

		return false;
	};

	std::cout << R"(
  [Space] Start/Pause simulation
  R,r     Toggle playback looping
  A,a     Toggle gravity
)";
	viewer.launch();
}
