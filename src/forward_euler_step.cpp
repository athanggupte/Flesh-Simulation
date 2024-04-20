#include "forward_euler_step.h"

namespace flesh {

	void forward_euler_step(
		Eigen::MatrixXd const& X,
		Eigen::MatrixXd const& V,
		Eigen::MatrixXd const& f,
		double h,
		Eigen::MatrixXd& X_next,
		Eigen::MatrixXd& V_next)
	{
		V_next = V + h * f;
		X_next = X + h * V_next;
	}

}
