#include "math_pseudoinverse_solver.h"

using namespace std;

namespace math {

// public func.

const char *PseudoinverseSolver::tag()
{
    return "PseudoinverseSolver";
}

PseudoinverseSolver::PseudoinverseSolver()
    :LinearSystemSolver()
{
}

PseudoinverseSolver::~PseudoinverseSolver()
{
}

std::string PseudoinverseSolver::id() const
{
    return std::string(PseudoinverseSolver::tag());
}

math::VectorNd_t PseudoinverseSolver::Solve(
        const math::MatrixN_t &coef_mat,
        const math::VectorNd_t &desired_vector
        ) const
{//TO DO	 
	/*
	Eigen::JacobiSVD<Eigen::MatrixXd> jacobian_svd(coef_mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
	return Eigen::VectorXd theta = jacobian_svd.solve(desired_vector);
	*/
	return math::VectorNd_t(coef_mat.transpose() * (coef_mat * coef_mat.transpose()).inverse() * desired_vector);
}

// protected func.

// private func.

} // namespace math {
