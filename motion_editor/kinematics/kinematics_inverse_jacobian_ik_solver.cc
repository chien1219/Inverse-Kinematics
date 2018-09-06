#include "kinematics_inverse_jacobian_ik_solver.h"
#include <limits>
#include "console_log.h"
#include "math_utils.h"
#include "math_linear_system_solver.h"
#include "acclaim_skeleton.h"
#include "acclaim_motion.h"
#include "kinematics_forward_solver.h"
#include "kinematics_pose.h"
#include "math_pseudoinverse_solver.h"

using namespace std;

namespace kinematics {

	InverseJacobianIkSolver::InverseJacobianIkSolver()
		:skeleton_(nullptr),
		fk_solver_(new ForwardSolver),
		step_(double{ 0.0 }),
		distance_epsilon_(double{ 0.0 }),
		max_iteration_num_(0),
		linear_system_solver_(nullptr)
	{
	}

	InverseJacobianIkSolver::~InverseJacobianIkSolver()
	{
	}

	void InverseJacobianIkSolver::Configure(
		const std::shared_ptr<acclaim::Skeleton> &skeleton,
		const std::shared_ptr<math::LinearSystemSolver> &linear_system_solver,
		const double step,
		const double distance_epsilon,
		const int32_t max_iteration_num
		)
	{
		skeleton_ = skeleton;
		fk_solver_->set_skeleton(skeleton_);
		fk_solver_->ConstructArticPath();

		linear_system_solver_ = linear_system_solver;

		step_ = step;
		distance_epsilon_ = distance_epsilon;
		max_iteration_num_ = max_iteration_num;
	}

	math::Vector6dColl_t InverseJacobianIkSolver::Solve(
		const math::Vector3d_t &target_pos,
		const int32_t start_bone_idx,
		const int32_t end_bone_idx,
		const math::Vector6dColl_t &original_whole_body_joint_pos6d
		)
	{
		//TO DO
		math::Vector6dColl_t new_pos6d = original_whole_body_joint_pos6d;
		PoseColl_t pose = fk_solver_->ComputeSkeletonPose(new_pos6d);
		int dof_cnt = 0;
		//Sum up dof
		int idx = start_bone_idx;
		while (idx != skeleton_->bone_ptr(end_bone_idx)->child->idx){
			dof_cnt += skeleton_->bone_ptr(idx)->dof;
			idx = skeleton_->bone_ptr(idx)->child->idx;
		}

		math::MatrixN_t Jac(3, dof_cnt);
		math::Vector3d_t pos = pose.at(end_bone_idx).end_pos();
		math::Vector3d_t dPos = (target_pos - pos) / max_iteration_num_;
		math::Vector3d_t tPos = pos;
		//Main Iteration Loop
		for (int i = 0; i < max_iteration_num_; i++){
			tPos += dPos;

			//Detect Job ends when close enough
			math::Vector3d_t V = target_pos - pos;
			if (sqrt(V.x() * V.x() + V.y() * V.y() + V.z() * V.z()) < distance_epsilon_)
				break;

			//Form Jacobian
			CalculateJacobian(tPos, start_bone_idx, end_bone_idx, Jac, pose);

			//Pseudo Inverse
			math::VectorNd_t newAngular = linear_system_solver_->Solve(Jac, tPos - pos);

			//Calculate and assign
			UpdateAngular(start_bone_idx, end_bone_idx, newAngular, new_pos6d);
			pose = fk_solver_->ComputeSkeletonPose(new_pos6d);
			pos = pose.at(end_bone_idx).end_pos();
		}
		return new_pos6d;
	}

	void InverseJacobianIkSolver::CalculateJacobian(
		const math::Vector3d_t &target_pos,
		const int32_t start_bone_idx,
		const int32_t end_bone_idx,
		math::MatrixN_t& Jac,
		PoseColl_t &pose){

		//Calculate Jacobian
		int idx = end_bone_idx;
		int i = 0;

		while (idx != skeleton_->bone_ptr(start_bone_idx)->parent->idx){

			math::Vector3d_t PR = target_pos - pose.at(idx).start_pos();
			//J = a X (p-r)
			math::Vector3d_t J;
			if (skeleton_->bone_ptr(idx)->dofx){
				math::Vector3d_t a = (pose.at(idx).rotation() * math::Vector3d_t(1, 0, 0)).normalized();
				Jac.col(i) = a.cross(PR);
				i++;
			}
			if (skeleton_->bone_ptr(idx)->dofy){
				math::Vector3d_t a = (pose.at(idx).rotation() * math::Vector3d_t(0, 1, 0)).normalized();
				Jac.col(i) = a.cross(PR);
				i++;
			}
			if (skeleton_->bone_ptr(idx)->dofz){
				math::Vector3d_t a = (pose.at(idx).rotation() * math::Vector3d_t(0, 0, 1)).normalized();
				Jac.col(i) = a.cross(PR);
				i++;
			}
			idx = skeleton_->bone_ptr(idx)->parent->idx;
		}
	}

	void InverseJacobianIkSolver::UpdateAngular(
		const int32_t start_bone_idx,
		const int32_t end_bone_idx,
		math::VectorNd_t &newAngular,
		math::Vector6dColl_t &new_pos6d){

		int i = 0;
		int idx = end_bone_idx;

		while (idx != skeleton_->bone_ptr(start_bone_idx)->parent->idx){

			math::Vector3d_t tmpAng = new_pos6d.at(idx).angular_vector();

			if (skeleton_->bone_ptr(idx)->dofx){
				tmpAng.x() += step_ * newAngular[i];
				i++;
			}
			if (skeleton_->bone_ptr(idx)->dofy){
				tmpAng.y() += step_ * newAngular[i];
				i++;
			}
			if (skeleton_->bone_ptr(idx)->dofz){
				tmpAng.z() += step_ * newAngular[i];
				i++;
			}
			new_pos6d.at(idx).set_angular_vector(tmpAng);
			idx = skeleton_->bone_ptr(idx)->parent->idx;
		}
	}

} // namespace kinematics {

