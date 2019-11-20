// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

#include <ceres/ceres.h>
#include "ceres/problem.h"
#include "ceres/solver.h"
#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/geometry/Similarity3_Kernel.hpp"
//- Robust estimation - LMeds (since no threshold can be defined)
#include "openMVG/robust_estimation/robust_estimator_LMeds.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

#include <ceres/rotation.h>
#include <ceres/types.h>

#include <iostream>
#include <limits>

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;

// Ceres CostFunctor used for SfM pose center to GPS pose center minimization
struct PoseCenterConstraintCostFunction
{
  Vec3 weight_;
  Vec3 pose_center_constraint_;

  PoseCenterConstraintCostFunction
  (
    const Vec3 & center,
    const Vec3 & weight
  ): weight_(weight), pose_center_constraint_(center)
  {
  }

  template <typename T> bool
  operator()
  (
    const T* const cam_extrinsics, // R_t
    T* residuals
  )
  const
  {
    using Vec3T = Eigen::Matrix<T,3,1>;
    Eigen::Map<const Vec3T> cam_R(&cam_extrinsics[0]);
    Eigen::Map<const Vec3T> cam_t(&cam_extrinsics[3]);
    const Vec3T cam_R_transpose(-cam_R);

    Vec3T pose_center;
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R_transpose.data(), cam_t.data(), pose_center.data());
    pose_center = pose_center * T(-1);

    Eigen::Map<Vec3T> residuals_eigen(residuals);
    residuals_eigen = weight_.cast<T>().cwiseProduct(pose_center - pose_center_constraint_.cast<T>());

    return true;
  }
};

struct ResidualErrorFunctor_Group_Intrinsic
{
	explicit ResidualErrorFunctor_Group_Intrinsic(const double* const pos_2dpoint)
		:m_pos_2dpoint(pos_2dpoint)
	{
	}

	// Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
	enum : uint8_t {
		OFFSET_FOCAL_LENGTH = 0,
		OFFSET_PRINCIPAL_POINT_X = 1,
		OFFSET_PRINCIPAL_POINT_Y = 2,
		OFFSET_DISTO_K1 = 3,
		OFFSET_DISTO_K2 = 4,
		OFFSET_DISTO_K3 = 5,
		OFFSET_DISTO_K4 = 6,
		//OFFSET_DISTO_T1 = 6,
		//OFFSET_DISTO_T2 = 7,
	};


	/**
	* @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
	* @param[in] cam_calibration:Camera calibration using one block of 6 parameters[R;t]
	* @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
	*   - 3 for rotation(angle axis), 3 for translation
	* @param[in] pos_3dpoint
	* @param[out] out_residuals
	*/
	template <typename T>
	bool operator()(
		const T* const cam_intrinsics,
		const T* const cam_calibration,
		const T* const cam_extrinsics,
		const T* const pos_3dpoint,
		T* out_residuals) const
	{
		//--
		// Apply external parameters (Pose)
		//--
		const T * cam_cR = cam_calibration;
		Eigen::Map<const Eigen::Matrix<T, 3, 1>> cam_ct(&cam_calibration[3]);
		Eigen::Matrix<T, 3, 3> R_i0;
		ceres::AngleAxisToRotationMatrix(cam_cR, R_i0.data());

		const T * cam_R = cam_extrinsics;
		Eigen::Map<const Eigen::Matrix<T, 3, 1>> cam_t(&cam_extrinsics[3]);
		Eigen::Matrix<T, 3, 3> R_0;
		ceres::AngleAxisToRotationMatrix(cam_R, R_0.data());

		Eigen::Matrix<T, 3, 3> R_i = R_i0*R_0;
		T angleAxis[3];
		ceres::RotationMatrixToAngleAxis(R_i.data(), angleAxis);

		Eigen::Matrix<T, 3, 1> transformed_point;
		// Rotate the point according the camera rotation
		ceres::AngleAxisRotatePoint(angleAxis, pos_3dpoint, transformed_point.data());

		// Apply the camera translation
		transformed_point += R_i*(R_0.transpose()*cam_t - cam_ct);

		// Transform the point from homogeneous to euclidean (undistorted point)
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		//--
		// Apply intrinsic parameters
		//--

		const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
		const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
		const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
		const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
		const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
		const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
		const T& k4 = cam_intrinsics[OFFSET_DISTO_K4];
		//const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
		//const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

		// Apply distortion (xd,yd) = disto(x_u,y_u)
		//鱼眼相机
		const T r2 = projected_point.squaredNorm();
		const T r = sqrt(r2);
		const T
			theta = atan(r),
			theta2 = theta*theta,
			theta3 = theta2*theta,
			theta4 = theta2*theta2,
			theta5 = theta4*theta,
			theta7 = theta3*theta3*theta, //thetha6*theta
			theta8 = theta4*theta4,
			theta9 = theta8*theta;
		const T theta_dist = theta + k1*theta3 + k2*theta5 + k3*theta7 + k4*theta9;
		const T inv_r = r > T(1e-8) ? T(1.0) / r : T(1.0);
		const T cdist = r > T(1e-8) ? theta_dist * inv_r : T(1.0);
		
		//布朗模型
		/*const T x_u = projected_point.x();
		const T y_u = projected_point.y();
		const T r2 = projected_point.squaredNorm();
		const T r4 = r2 * r2;
		const T r6 = r4 * r2;
		const T r_coeff = (1.0 + k1 * r2 + k2 * r4 + k3 * r6);
		const T t_x = t2 * (r2 + 2.0 * x_u * x_u) + 2.0 * t1 * x_u * y_u;
		const T t_y = t1 * (r2 + 2.0 * y_u * y_u) + 2.0 * t2 * x_u * y_u;*/

		// Compute and return the error is the difference between the predicted
		//  and observed position
		Eigen::Map<Eigen::Matrix<T, 2, 1>> residuals(out_residuals);
		//鱼眼
		residuals << principal_point_x + (projected_point.x() * cdist) * focal - m_pos_2dpoint[0],
			principal_point_y + (projected_point.y() * cdist) * focal - m_pos_2dpoint[1];

		//布朗模型
		/*residuals << principal_point_x + (projected_point.x() * r_coeff + t_x) * focal - m_pos_2dpoint[0],
			principal_point_y + (projected_point.y() * r_coeff + t_y) * focal - m_pos_2dpoint[1];*/
		return true;
	}

	static int num_residuals() { return 2; }

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create
	(
		const Vec2 & observation,
		const double weight = 0.0
	)
	{
		if (weight == 0.0)
		{
			return
				(new ceres::AutoDiffCostFunction
					<ResidualErrorFunctor_Group_Intrinsic, 2,7, 6, 6, 3>(
						new ResidualErrorFunctor_Group_Intrinsic(observation.data())));
		}
		else
		{
			return
				(new ceres::AutoDiffCostFunction
					<WeightedCostFunction_group<ResidualErrorFunctor_Group_Intrinsic>, 2, 7, 6, 6, 3>
					(new WeightedCostFunction_group<ResidualErrorFunctor_Group_Intrinsic>
					(new ResidualErrorFunctor_Group_Intrinsic(observation.data()), weight)));
		}
	}

	const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Group_Intrinsic_Multi
{
	explicit ResidualErrorFunctor_Group_Intrinsic_Multi(const double* const pos_2dpoint)
		:m_pos_2dpoint(pos_2dpoint)
	{
	}

	// Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
	enum : uint8_t {
		OFFSET_FOCAL_LENGTH = 0,
		OFFSET_PRINCIPAL_POINT_X = 1,
		OFFSET_PRINCIPAL_POINT_Y = 2,
		OFFSET_DISTO_K1 = 3,
		OFFSET_DISTO_K2 = 4,
		OFFSET_DISTO_K3 = 5,
		OFFSET_DISTO_K4 = 6,
		//OFFSET_DISTO_T1 = 6,
		//OFFSET_DISTO_T2 = 7,
	};

	template <typename T>
	bool operator()(
		const T* const cam_intrinsics,
		const T* const cam_calibration,
		const T* const cam_sevenparameters,
		const T* const cam_extrinsics,
		const T* const pos_3dpoint,
		T* out_residuals) const
	{
		//--
		// Apply external parameters (Pose)
		//--
		const T * cam_cR = cam_calibration;
		Eigen::Map<const Eigen::Matrix<T, 3, 1>> cam_ct(&cam_calibration[3]);
		Eigen::Matrix<T, 3, 3> R_i0;
		ceres::AngleAxisToRotationMatrix(cam_cR, R_i0.data());

		const T * cam_cS = cam_sevenparameters;
		Eigen::Map<const Eigen::Matrix<T, 3, 1>> cam_cSt(&cam_sevenparameters[3]);
		Eigen::Matrix<T, 3, 3> R_S;
		ceres::AngleAxisToRotationMatrix(cam_cS, R_S.data());
		//将pos_3dpoint的坐标系转换
		Eigen::Map<const Eigen::Matrix<T, 3, 1>> pos_3d_point(&pos_3dpoint[0]);
		Eigen::Matrix<T, 3, 1> pos_3dpointS = R_S*(pos_3d_point - cam_cSt);

		const T * cam_R = cam_extrinsics;
		Eigen::Map<const Eigen::Matrix<T, 3, 1>> cam_t(&cam_extrinsics[3]);
		Eigen::Matrix<T, 3, 3> R_0;
		ceres::AngleAxisToRotationMatrix(cam_R, R_0.data());

		Eigen::Matrix<T, 3, 3> R_i = R_i0*R_0;
		T angleAxis[3];
		ceres::RotationMatrixToAngleAxis(R_i.data(), angleAxis);

		Eigen::Matrix<T, 3, 1> transformed_point;
		// Rotate the point according the camera rotation
		ceres::AngleAxisRotatePoint(angleAxis, pos_3dpointS.data(), transformed_point.data());

		// Apply the camera translation
		transformed_point += R_i*(R_0.transpose()*cam_t - cam_ct);

		// Transform the point from homogeneous to euclidean (undistorted point)
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		//--
		// Apply intrinsic parameters
		//--

		const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
		const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
		const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
		const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
		const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
		const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
		const T& k4 = cam_intrinsics[OFFSET_DISTO_K4];
		//const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
		//const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

		// Apply distortion (xd,yd) = disto(x_u,y_u)
		//鱼眼相机
		const T r2 = projected_point.squaredNorm();
		const T r = sqrt(r2);
		const T
			theta = atan(r),
			theta2 = theta*theta,
			theta3 = theta2*theta,
			theta4 = theta2*theta2,
			theta5 = theta4*theta,
			theta7 = theta3*theta3*theta, //thetha6*theta
			theta8 = theta4*theta4,
			theta9 = theta8*theta;
		const T theta_dist = theta + k1*theta3 + k2*theta5 + k3*theta7 + k4*theta9;
		const T inv_r = r > T(1e-8) ? T(1.0) / r : T(1.0);
		const T cdist = r > T(1e-8) ? theta_dist * inv_r : T(1.0);

		//布朗模型
		/*const T x_u = projected_point.x();
		const T y_u = projected_point.y();
		const T r2 = projected_point.squaredNorm();
		const T r4 = r2 * r2;
		const T r6 = r4 * r2;
		const T r_coeff = (1.0 + k1 * r2 + k2 * r4 + k3 * r6);
		const T t_x = t2 * (r2 + 2.0 * x_u * x_u) + 2.0 * t1 * x_u * y_u;
		const T t_y = t1 * (r2 + 2.0 * y_u * y_u) + 2.0 * t2 * x_u * y_u;*/

		// Compute and return the error is the difference between the predicted
		//  and observed position
		Eigen::Map<Eigen::Matrix<T, 2, 1>> residuals(out_residuals);
		//鱼眼
		residuals << principal_point_x + (projected_point.x() * cdist) * focal - m_pos_2dpoint[0],
			principal_point_y + (projected_point.y() * cdist) * focal - m_pos_2dpoint[1];

		//布朗模型
		/*residuals << principal_point_x + (projected_point.x() * r_coeff + t_x) * focal - m_pos_2dpoint[0],
		principal_point_y + (projected_point.y() * r_coeff + t_y) * focal - m_pos_2dpoint[1];*/
		return true;
	}

	static int num_residuals() { return 2; }

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create
	(
		const Vec2 & observation,
		const double weight = 0.0
	)
	{
		if (weight == 0.0)
		{
			return
				(new ceres::AutoDiffCostFunction
					<ResidualErrorFunctor_Group_Intrinsic_Multi, 2, 7, 6, 6, 6, 3>(
						new ResidualErrorFunctor_Group_Intrinsic_Multi(observation.data())));
		}
		else
		{
			return
				(new ceres::AutoDiffCostFunction
					<WeightedCostFunction_group<ResidualErrorFunctor_Group_Intrinsic_Multi>, 2, 7, 6, 6, 6, 3>
					(new WeightedCostFunction_group<ResidualErrorFunctor_Group_Intrinsic_Multi>
					(new ResidualErrorFunctor_Group_Intrinsic_Multi(observation.data()), weight)));
		}
	}

	const double * m_pos_2dpoint; // The 2D observation
};
/// Create the appropriate cost functor according the provided input camera intrinsic model.
/// The residual can be weighetd if desired (default 0.0 means no weight).
ceres::CostFunction * IntrinsicsToCostFunction
(
  IntrinsicBase * intrinsic,
  const Vec2 & observation,
  const double weight
)
{
  switch (intrinsic->getType())
  {
    case PINHOLE_CAMERA:
      return ResidualErrorFunctor_Pinhole_Intrinsic::Create(observation, weight);
    case PINHOLE_CAMERA_RADIAL1:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1::Create(observation, weight);
    case PINHOLE_CAMERA_RADIAL3:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3::Create(observation, weight);
    case PINHOLE_CAMERA_BROWN:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2::Create(observation, weight);
    case PINHOLE_CAMERA_FISHEYE:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye::Create(observation, weight);
    case CAMERA_SPHERICAL:
      return ResidualErrorFunctor_Intrinsic_Spherical::Create(intrinsic, observation, weight);
    default:
      return {};
  }
}

ceres::CostFunction * MaincamIntrinsicsToCostFunction
(
	IntrinsicBase * intrinsic,
	const Vec2 & observation,
	const double weight
)
{
	//return ResidualErrorFunctor_Group_Intrinsic::Create(observation, weight);
	return ResidualErrorFunctor_Group_Intrinsic_Multi::Create(observation, weight);
}

Bundle_Adjustment_Ceres::BA_Ceres_options::BA_Ceres_options
(
  const bool bVerbose,
  bool bmultithreaded
)
: bVerbose_(bVerbose),
  nb_threads_(1),
  parameter_tolerance_(1e-8), //~= numeric_limits<float>::epsilon()
  bUse_loss_function_(true)
{
  #ifdef OPENMVG_USE_OPENMP
    nb_threads_ = omp_get_max_threads();
  #endif // OPENMVG_USE_OPENMP
  if (!bmultithreaded)
    nb_threads_ = 1;

  bCeres_summary_ = false;

  // Default configuration use a DENSE representation
  linear_solver_type_ = ceres::DENSE_SCHUR;
  preconditioner_type_ = ceres::JACOBI;
  // If Sparse linear solver are available
  // Descending priority order by efficiency (SUITE_SPARSE > CX_SPARSE > EIGEN_SPARSE)
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
  {
    sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
    linear_solver_type_ = ceres::SPARSE_SCHUR;
  }
  else
  {
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
    {
      sparse_linear_algebra_library_type_ = ceres::CX_SPARSE;
      linear_solver_type_ = ceres::SPARSE_SCHUR;
    }
    else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
    {
      sparse_linear_algebra_library_type_ = ceres::EIGEN_SPARSE;
      linear_solver_type_ = ceres::SPARSE_SCHUR;
    }
  }
}


Bundle_Adjustment_Ceres::Bundle_Adjustment_Ceres
(
  const Bundle_Adjustment_Ceres::BA_Ceres_options & options
)
: ceres_options_(options)
{}

Bundle_Adjustment_Ceres::BA_Ceres_options &
Bundle_Adjustment_Ceres::ceres_options()
{
  return ceres_options_;
}

bool Bundle_Adjustment_Ceres::Adjust
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------

  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  std::vector<std::pair<int,Vec3>> Weight_XYZ;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
	  std::cout << "options.use_motion_priors_opt=" << options.use_motion_priors_opt << std::endl;
	// Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
		  Weight_XYZ.push_back(std::make_pair(prior->id_pose, prior->center_weight_));
        }
      }

	  if(options.use_sim_opt)
	  	{
	  	  openMVG::geometry::Similarity3 sim;
		  // Compute the registration:
		  if (X_GPS.size() > 3)
		  {
			  const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(), 3, X_SfM.size());
			  const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(), 3, X_GPS.size());
			  geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
			  const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);

			  if (lmeds_median != std::numeric_limits<double>::max())
			  {
				  std::cerr << "Compute the registration:lmeds_median(最小中值二乘)=" << lmeds_median << std::endl;
				  b_usable_prior = true; // PRIOR can be used safely

				 // Compute the median residual error once the registration is applied
				  for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
				  {
					  pos = sim(pos);
				  }
			  }
		  }
		  // Apply the found transformation to the SfM Data Scene
		openMVG::sfm::ApplySimilarity(sim, sfm_data);
	  }

	 if (X_GPS.size() > 3)
	  {
		  b_usable_prior = true;
		  Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
		  std::sort(residual.data(), residual.data() + residual.size());
		  //残差中位数
		  pose_center_robust_fitting_error = residual(residual.size() / 2);
		  ////均方根误差
		  double average_fitting_error = 0;
		  double sum_error = 0;
		  double square_error = 0;
		  for (int i = 0; i < residual.size(); i++)
			  sum_error += residual[i];
		  average_fitting_error = sum_error / residual.size();
		  sum_error = 0;
		  for (int i = 0; i < residual.size(); i++)
			  sum_error += std::pow(residual[i] - average_fitting_error, 2);
		  const double root_mean_square_error = sqrt(sum_error / residual.size());
		  
		  //根据residual值和lmeds_median值加权
		  for (int i = 0; i < residual.size(); i++)
		  {
			  if (residual[i] >= root_mean_square_error * 3)
				  Weight_XYZ[i].second = Vec3(0, 0, 0);
			  if (residual[i] <= root_mean_square_error)
				  Weight_XYZ[i].second = Vec3(1, 1, 1);
			  if (residual[i] > root_mean_square_error&residual[i] < 3 * root_mean_square_error)
				  Weight_XYZ[i].second = Vec3(pow(residual[i] / (3 * root_mean_square_error), 2), pow(residual[i] / (3 * root_mean_square_error), 2), pow(residual[i] / (3 * root_mean_square_error), 2));
		  }

		  // Move entire scene to center for better numerical stability
		  Vec3 pose_centroid = Vec3::Zero();
		  for (const auto & pose_it : sfm_data.poses)
		  {
			  pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
		  }
		  sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
		  openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
	  }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double>> map_intrinsics;
  Hash_Map<IndexT, std::vector<double>> map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses.at(indexPose)[0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.insert(vec_constant_extrinsic.end(), {0,1,2});
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.insert(vec_constant_extrinsic.end(), {3,4,5});
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics.at(indexCam).empty())
      {
        double * parameter_block = &map_intrinsics.at(indexCam)[0];
        problem.AddParameterBlock(parameter_block, map_intrinsics.at(indexCam).size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics.at(indexCam).size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;
  
  //Use Control Point
  //std::set<int> cp_ids = { 0,1,3,4,5,6,10 };
  if (options.control_point_opt.bUse_control_points)
  {
	  // Use Ground Control Point:
	  
	  //坐标转换
	  {
		  std::map<IndexT, Vec3> vec_control_points, vec_triangulated;
		  std::map<IndexT, double> vec_triangulation_errors;
		  for (const auto & control_point_it : sfm_data.control_points)
		  {
			  const Landmark & landmark = control_point_it.second;
			  //Triangulate the observations:
			  const Observations & obs = landmark.obs;
			  std::vector<Vec3> bearing;
			  std::vector<Mat34> poses;
			  bearing.reserve(obs.size());
			  poses.reserve(obs.size());
			  for (const auto & obs_it : obs)
			  {
				  const View * view = sfm_data.views.at(obs_it.first).get();
				  if (!sfm_data.IsPoseAndIntrinsicDefined(view))
					  continue;
				  const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
				  const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
				  const Vec2 pt = obs_it.second.x;
				  bearing.emplace_back((*cam)(cam->get_ud_pixel(pt)));
				  poses.emplace_back(pose.asMatrix());
			  }
			  const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
			  Vec4 Xhomogeneous;
			  TriangulateNView(bearing_matrix, poses, &Xhomogeneous);
			  const Vec3 X = Xhomogeneous.hnormalized();
			  // Test validity of the hypothesis (front of the cameras):
			  bool bChierality = true;
			  int i(0);
			  double reprojection_error_sum(0.0);
			  for (const auto & obs_it : obs)
			  {
				  const View * view = sfm_data.views.at(obs_it.first).get();
				  if (!sfm_data.IsPoseAndIntrinsicDefined(view))
					  continue;

				  const Pose3 pose = sfm_data.GetPoseOrDie(view);
				  bChierality &= CheiralityTest(bearing[i], pose, X);
				  const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
				  const Vec2 pt = obs_it.second.x;
				  const Vec2 residual = cam->residual(pose(X), pt);
				  reprojection_error_sum += residual.norm();
				  ++i;
			  }
			  if (bChierality) // Keep the point only if it has a positive depth
			  {
				  vec_triangulated[control_point_it.first] = X;
				  vec_control_points[control_point_it.first] = landmark.X;
				  vec_triangulation_errors[control_point_it.first] = reprojection_error_sum / (double)bearing.size();
			  }
			  else
			  {
				  std::cout << "Control Point cannot be triangulated (not in front of the cameras)" << std::endl;
				  return false;
			  }
		  }

		  if (vec_control_points.size() < 3)
		  {
			  std::cout << "Insufficient number of triangulated control points." << std::endl;
			  return false;
		  }

		  // compute the similarity
		  {
			  // data conversion to appropriate container
			  Mat x1(3, vec_control_points.size()),
				  x2(3, vec_control_points.size());
			  for (size_t i = 0; i < vec_control_points.size(); ++i)
			  {
				  x1.col(i) = vec_triangulated[i];
				  x2.col(i) = vec_control_points[i];
			  }

			  std::cout
				  << "Control points observation triangulations:\n"
				  << x1 << std::endl << std::endl
				  << "Control points coords:\n"
				  << x2 << std::endl << std::endl;

			  Vec3 t;
			  Mat3 R;
			  double S;
			  if (openMVG::geometry::FindRTS(x1, x2, &S, &t, &R))
			  {
				  openMVG::geometry::Refine_RTS(x1, x2, &S, &t, &R);
				  std::cout << "Found transform:\n"
					  << " scale: " << S << "\n"
					  << " rotation:\n" << R << "\n"
					  << " translation: " << t.transpose() << std::endl;
				  //--
				  // Apply the found transformation as a 3D Similarity transformation matrix // S * R * X + t
				  //--

				  const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t / S), S);
				  openMVG::sfm::ApplySimilarity(sim, sfm_data);
			  }
			  else
			  {
				  std::cout << "Registration failed. Please check your Control Points coordinates." << std::endl;
				  return false;
			  }
			  vec_control_points.clear();
			  vec_triangulated.clear();
			  vec_triangulation_errors.clear();
		  }
	  }
	  
	  for (auto & gcp_landmark_it : sfm_data.control_points)
	  {
		  //if(cp_ids.count(gcp_landmark_it.first)==0)
			 // continue;
		  const Observations & obs = gcp_landmark_it.second.obs;

		  for (const auto & obs_it : obs)
		  {
			  // Build the residual block corresponding to the track observation:
			  const View * view = sfm_data.views.at(obs_it.first).get();

			  // Each Residual block takes a point and a camera as input and outputs a 2
			  // dimensional residual. Internally, the cost function stores the observed
			  // image location and compares the reprojection against the observation.
			  ceres::CostFunction* cost_function =
				  IntrinsicsToCostFunction(
					  sfm_data.intrinsics.at(view->id_intrinsic).get(),
					  obs_it.second.x,
					  options.control_point_opt.weight);

			  if (cost_function)
			  {
				  if (!map_intrinsics.at(view->id_intrinsic).empty())
				  {
					  problem.AddResidualBlock(cost_function,
						  nullptr,
						  &map_intrinsics.at(view->id_intrinsic)[0],
						  &map_poses.at(view->id_pose)[0],
						  gcp_landmark_it.second.X.data());
				  }
				  else
				  {
					  problem.AddResidualBlock(cost_function,
						  nullptr,
						  &map_poses.at(view->id_pose)[0],
						  gcp_landmark_it.second.X.data());
				  }
			  }
		  }
		  if (obs.empty())
		  {
			  std::cerr
				  << "Cannot use this GCP id: " << gcp_landmark_it.first
				  << ". There is not linked image observation." << std::endl;
		  }
		  else
		  {
			  // Set the 3D point as FIXED (it's a valid GCP)
			  problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
		  }
	  }
  }

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
	  const Observations & obs = structure_landmark_it.second.obs;

	  for (const auto & obs_it : obs)
	  {
		  // Build the residual block corresponding to the track observation:
		  const View * view = sfm_data.views.at(obs_it.first).get();

		  // Each Residual block takes a point and a camera as input and outputs a 2
		  // dimensional residual. Internally, the cost function stores the observed
		  // image location and compares the reprojection against the observation.
		  ceres::CostFunction* cost_function =
			  IntrinsicsToCostFunction(sfm_data.intrinsics.at(view->id_intrinsic).get(),
				  obs_it.second.x);

		  if (cost_function)
		  {
			  if (!map_intrinsics.at(view->id_intrinsic).empty())
			  {
				  problem.AddResidualBlock(cost_function,
					  p_LossFunction,
					  &map_intrinsics.at(view->id_intrinsic)[0],
					  &map_poses.at(view->id_pose)[0],
					  structure_landmark_it.second.X.data());
			  }
			  else
			  {
				  problem.AddResidualBlock(cost_function,
					  p_LossFunction,
					  &map_poses.at(view->id_pose)[0],
					  structure_landmark_it.second.X.data());
			  }
		  }
		  else
		  {
			  std::cerr << "Cannot create a CostFunction for this camera model." << std::endl;
			  return false;
		  }
	  }
	  if (options.structure_opt == Structure_Parameter_Type::NONE)
		  problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
		  Vec3 weight;
		  //认为XYZ的权值一样
		  for (const auto & it_w:Weight_XYZ)
			   if (it_w.first == prior->id_pose)
				   weight = it_w.second[0]* prior->center_weight_;
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_,weight));

        problem.AddResidualBlock(
          cost_function,
          new ceres::HuberLoss(
            Square(pose_center_robust_fitting_error)),
                   &map_poses.at(prior->id_view)[0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
#if CERES_VERSION_MAJOR < 2
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
#endif
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses.at(indexPose)[0], R_refined.data());
        Vec3 t_refined(map_poses.at(indexPose)[3], map_poses.at(indexPose)[4], map_poses.at(indexPose)[5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics.at(indexCam);
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    //检查精度
	if (options.control_point_opt.bUse_control_points)
	{
		for (auto & cp_it : sfm_data.control_points)
		{
			int index_cpt = cp_it.first;
			openMVG::sfm::Landmark& check_point = sfm_data.control_points[index_cpt];
			std::cout << "the coordinate of the check point(" << index_cpt << "):" << check_point.X[0] << " " << check_point.X[1] << " " << check_point.X[2] << std::endl;
			//Triangulate the observations:
			const openMVG::sfm::Observations & cobs = check_point.obs;
			std::cout << "The number of check_point's visibility:" << cobs.size() << std::endl;
			std::vector<openMVG::Vec3> cbearing;
			std::vector<openMVG::Mat34> cposes;
			double reprojection_error_sum(0);
			cbearing.reserve(cobs.size());
			cposes.reserve(cobs.size());
			for (const auto & cobs_it : cobs)
			{
				const openMVG::sfm::View * view = sfm_data.views.at(cobs_it.first).get();
				std::cout << "View ID:" << view->id_view << "  ";
				if (!sfm_data.IsPoseAndIntrinsicDefined(view))
					continue;
				const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
				const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
				const openMVG::Vec2 pt = cobs_it.second.x;
				std::cout << pt[0] << " " << pt[1];
				cbearing.emplace_back((*cam)(cam->get_ud_pixel(pt)));
				cposes.emplace_back(pose.asMatrix());
				std::cout << std::endl;
			}
			const Eigen::Map<const openMVG::Mat3X> cbearing_matrix(cbearing[0].data(), 3, cbearing.size());
			openMVG::Vec4 cXhomogeneous;

			openMVG::TriangulateNView(cbearing_matrix, cposes, &cXhomogeneous);

			const openMVG::Vec3 cX = cXhomogeneous.hnormalized();

			//三角化的重投影误差
			int index_coord = 0;
			for (const auto & cobs_it : cobs)
			{
				index_coord = cobs_it.second.id_feat;
				const openMVG::sfm::View * view = sfm_data.views.at(cobs_it.first).get();
				if (!sfm_data.IsPoseAndIntrinsicDefined(view))
					continue;
				const openMVG::cameras::IntrinsicBase * cam_ = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
				const openMVG::geometry::Pose3 pose_ = sfm_data.GetPoseOrDie(view);
				const openMVG::Vec2 pt_ = cobs_it.second.x;
				const Vec2 residual = cam_->residual(pose_(cX), pt_);
				reprojection_error_sum += residual.norm();
			}
			double mean_reprojection_error = reprojection_error_sum / (double)cobs.size();
			
			//对于控制点要进行坐标系转换
			openMVG::Vec3 dX = cX - check_point.X;
			std::cout << "the triangulated coordinate of the check point:" << cX[0] << " " << cX[1] << " " << cX[2] << std::endl;
			std::cout << "the mean reprojection error after BA:" << mean_reprojection_error << std::endl;
			std::cout << "the bias after BA:" << dX[0] << " " << dX[1] << " " << dX[2] << std::endl;
		}
	}

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_group
(
	SfM_Data & sfm_data,     // the SfM scene to refine
	std::vector<Pose3>& calibrations,
	const Optimize_Options & options
)
{
	//----------
	// Add camera parameters
	// - intrinsics
	// - poses [R|t]

	// Create residuals for each observation in the bundle adjustment problem. The
	// parameters for cameras and points are added automatically.
	//----------

	double pose_center_robust_fitting_error = 0.0;
	openMVG::geometry::Similarity3 sim_to_center;
	bool b_usable_prior = false;
	std::vector<std::pair<int, Vec3>> Weight_XYZ;
	std::map<int, geometry::Similarity3> seven_parameters;
	int reference_coors = 0;
	seven_parameters[reference_coors] = openMVG::geometry::Similarity3();
	std::cout <<"sfm_data.poses.size()"<< sfm_data.poses.size()<< ",calibrations.size()=" << calibrations.size() << std::endl;
	int station_num = sfm_data.poses.size() / calibrations.size();
	std::cout << "station_num=" << station_num << std::endl;
	if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
	{
		// - Compute a robust X-Y affine transformation & apply it
		// - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)

		// Collect corresponding camera centers
		std::vector<Vec3> X_SfM, X_GPS;
		for (const auto & view_it : sfm_data.GetViews())
		{
			const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
			if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
			{
				X_SfM.push_back(sfm_data.GetPoses().at(prior->id_pose).center());
				X_GPS.push_back(prior->pose_center_);
				Weight_XYZ.push_back(std::make_pair(prior->id_pose, prior->center_weight_));
			}
		}

		if (options.use_sim_opt)
		{
			openMVG::geometry::Similarity3 sim;
			// Compute the registration:
			if (X_GPS.size() > 3)
			{
				const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(), 3, X_SfM.size());
				const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(), 3, X_GPS.size());
				geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
				const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);

				if (lmeds_median != std::numeric_limits<double>::max())
				{
					std::cerr << "Compute the registration:lmeds_median(最小中值二乘)=" << lmeds_median << std::endl;
					b_usable_prior = true; // PRIOR can be used safely

										   // Compute the median residual error once the registration is applied
					for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
					{
						pos = sim(pos);
					}
				}
			}
			// Apply the found transformation to the SfM Data Scene
			openMVG::sfm::ApplySimilarity(sim, sfm_data);
		}

		if (X_GPS.size() > 3)
		{
			b_usable_prior = true;
			Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
			std::sort(residual.data(), residual.data() + residual.size());
			//残差中位数
			pose_center_robust_fitting_error = residual(residual.size() / 2);
			////均方根误差
			double average_fitting_error = 0;
			double sum_error = 0;
			double square_error = 0;
			for (int i = 0; i < residual.size(); i++)
				sum_error += residual[i];
			average_fitting_error = sum_error / residual.size();
			sum_error = 0;
			for (int i = 0; i < residual.size(); i++)
				sum_error += std::pow(residual[i] - average_fitting_error, 2);
			const double root_mean_square_error = sqrt(sum_error / residual.size());

			//根据residual值和lmeds_median值加权
			for (int i = 0; i < residual.size(); i++)
			{
				if (residual[i] >= root_mean_square_error * 3)
					Weight_XYZ[i].second = Vec3(0, 0, 0);
				if (residual[i] <= root_mean_square_error)
					Weight_XYZ[i].second = Vec3(1, 1, 1);
				if (residual[i] > root_mean_square_error&residual[i] < 3 * root_mean_square_error)
					Weight_XYZ[i].second = Vec3(pow(residual[i] / (3 * root_mean_square_error), 2), pow(residual[i] / (3 * root_mean_square_error), 2), pow(residual[i] / (3 * root_mean_square_error), 2));
			}

			// Move entire scene to center for better numerical stability
			Vec3 pose_centroid = Vec3::Zero();
			for (const auto & pose_it : sfm_data.poses)
			{
				pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
			}
			sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
			openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
		}
	}

	//控制点为多坐标系
	if (options.control_point_opt.bUse_control_points && options.control_point_opt.bUse_multi_coordinate)
	{
		//首先遍历所有的控制点判断有多少个坐标系
		//每个坐标系对应的控制点存到一个set里面
		std::map<int,std::set<int>>kzd_coords;

		for (auto & gcp_landmark_it : sfm_data.control_points)
		{
			const Observations & obs = gcp_landmark_it.second.obs;
			kzd_coords[obs.begin()->second.id_feat].insert(gcp_landmark_it.first);
		}
		
		//统计每个坐标系的观测值数
		std::map<int, int>coords_obs;
		for (auto coors_it=kzd_coords.begin();coors_it!=kzd_coords.end();++coors_it)
		{
			//删除控制点数目少于3的坐标系
			if (coors_it->second.size() < 3)
			{
				kzd_coords.erase(coors_it);
				continue;
			}
			int count = 0;
			//删除观测值少于2的控制点
			for (auto landmarkid_it : coors_it->second)
			{
				int obs_size= sfm_data.control_points.at(landmarkid_it).obs.size();
				if (obs_size < 2)
				{
					coors_it->second.erase(landmarkid_it);
					continue;
				}
				count += obs_size;
			}
				
			coords_obs[coors_it->first] = count;
			std::cout << "The coords_" << coors_it->first << " ,obs:" << coords_obs[coors_it->first] << ",control_points_num:" <<coors_it->second.size() << std::endl;
		}
			
		//然后以控制点平均观测值最多的一个坐标系为参考
		double avg_obs_max = 0.0;
		for (auto & refer_it : coords_obs)
		{
			double avg_obs = 1.0*refer_it.second / kzd_coords[refer_it.first].size();
			std::cout << "The coords_" << refer_it.first << " has average obs:" << avg_obs << std::endl;
			if (avg_obs > avg_obs_max)
			{
				reference_coors = refer_it.first;
				avg_obs_max = avg_obs;
			}
		}
		std::cout << "The reference coordinate is:" << reference_coors << std::endl;

		//参考坐标系的七参数
		seven_parameters[reference_coors] = geometry::Similarity3();
		
		//进行七参转换
		{
			std::cout << "..............................The first transforming..............................." << std::endl;
			std::map<IndexT, Vec3> vec_control_points, vec_triangulated;
			std::map<IndexT, double> vec_triangulation_errors;
			for (auto & landmarkid_it : kzd_coords[reference_coors])
			{
				std::cout << "landmarkid_it=" << landmarkid_it << std::endl;
				const openMVG::sfm::Landmark& landmark = sfm_data.control_points[landmarkid_it];
				const Observations & obs = landmark.obs;
				std::cout << "Triangulate the observations("<<obs.size()<<").................." << std::endl;
				std::vector<Vec3> bearing;
				std::vector<Mat34> poses;
				bearing.reserve(obs.size());
				poses.reserve(obs.size());
				for (const auto & obs_it : obs)
				{
					const View * view = sfm_data.views.at(obs_it.first).get();
					if (!sfm_data.IsPoseAndIntrinsicDefined(view))
						continue;
					const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
					const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
					const Vec2 pt = obs_it.second.x;
					bearing.emplace_back((*cam)(cam->get_ud_pixel(pt)));
					poses.emplace_back(pose.asMatrix());
				}
				const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
				Vec4 Xhomogeneous;
				TriangulateNView(bearing_matrix, poses, &Xhomogeneous);
				const Vec3 X = Xhomogeneous.hnormalized();
				// Test validity of the hypothesis (front of the cameras):
				bool bChierality = true;
				int i(0);
				double reprojection_error_sum(0.0);
				for (const auto & obs_it : obs)
				{
					const View * view = sfm_data.views.at(obs_it.first).get();
					if (!sfm_data.IsPoseAndIntrinsicDefined(view))
						continue;

					const Pose3 pose = sfm_data.GetPoseOrDie(view);
					bChierality &= CheiralityTest(bearing[i], pose, X);
					const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
					const Vec2 pt = obs_it.second.x;
					const Vec2 residual = cam->residual(pose(X), pt);
					reprojection_error_sum += residual.norm();
					++i;
				}
				if (bChierality) // Keep the point only if it has a positive depth
				{
					vec_triangulated[landmarkid_it] = X;
					std::cout << "vec_triangulated[" << landmarkid_it << "] =" << X[0]<<","<<X[1]<<","<<X[2]<< std::endl;
					vec_control_points[landmarkid_it] = landmark.X;
					std::cout << "vec_control_points[" << landmarkid_it << "] =" << landmark.X[0] << "," << landmark.X[1] << "," << landmark.X[2] << std::endl;
					vec_triangulation_errors[landmarkid_it] = reprojection_error_sum / (double)bearing.size();
				}
				else
				{
					std::cout << "Control Point cannot be triangulated (not in front of the cameras)" << std::endl;
					return false;
				}
			}
			if (vec_control_points.size() < 3)
			{
				std::cout << "Insufficient number of triangulated control points." << std::endl;
				return false;
			}

			// data conversion to appropriate container
			Mat x1(3, vec_control_points.size()),
				x2(3, vec_control_points.size());
			int col_i(0);
			for (auto & cp_it : vec_control_points)
			{
				x1.col(col_i) = vec_triangulated[cp_it.first];
				x2.col(col_i) = vec_control_points[cp_it.first];
				++col_i;
			}

			std::cout
				<< "Control points observation triangulations:\n"
				<< x1 << std::endl << std::endl
				<< "Control points coords:\n"
				<< x2 << std::endl << std::endl;
			Vec3 t;
			Mat3 R;
			double S;
			if (openMVG::geometry::FindRTS(x1, x2, &S, &t, &R))
			{
				openMVG::geometry::Refine_RTS(x1, x2, &S, &t, &R);
				std::cout << "Found transform:\n"
					<< " scale: " << S << "\n"
					<< " rotation:\n" << R << "\n"
					<< " translation: " << t.transpose() << std::endl;

				//--
				// Apply the found transformation as a 3D Similarity transformation matrix // S * R * X + t
				//--
				
				const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t / S), S);
				openMVG::sfm::ApplySimilarity(sim, sfm_data);

				// Display some statistics:
				std::stringstream os;
				for (auto & landmarkid_it : kzd_coords[reference_coors])
				{
					const IndexT CPIndex = landmarkid_it;
					// If the control point has not been used, continue...
					if (vec_triangulation_errors.find(CPIndex) == vec_triangulation_errors.end())
						continue;

					os
						<< "CP index: " << CPIndex << "\n"
						<< "CP triangulation error: " << vec_triangulation_errors[CPIndex] << " pixel(s)\n"
						<< "CP registration error: "
						<< (sim(vec_triangulated[CPIndex]) - vec_control_points[CPIndex]).norm() << " user unit(s)" << "\n\n";
				}
				std::cout << os.str();
				
			}
		}

		//计算其他坐标系的seven_parameters
		{
			for (auto & coors_it : kzd_coords)
			{
				std::cout << "!!!!!!!!!!!!!!Beginning to evaluate the seven parameters of Coordinate_" << coors_it.first << std::endl;
				if(coors_it.first==reference_coors)
					continue;
				//控制点数少于三的坐标系不参与计算
				if (coors_it.second.size() < 3)
				{
					std::cout << "Coordinate_id=" << coors_it.first << ",not enough control point(<3),deserted!";
					continue;
				}
				//七参计算
				{
					std::map<IndexT, Vec3> vec_control_points, vec_triangulated;
					std::map<IndexT, double> vec_triangulation_errors;
					for (auto & landmarkid_it : kzd_coords[coors_it.first])
					{
						const openMVG::sfm::Landmark& landmark = sfm_data.control_points[landmarkid_it];
						//Triangulate the observations:
						const Observations & obs = landmark.obs;
						std::vector<Vec3> bearing;
						std::vector<Mat34> poses;
						bearing.reserve(obs.size());
						poses.reserve(obs.size());
						for (const auto & obs_it : obs)
						{
							const View * view = sfm_data.views.at(obs_it.first).get();
							if (!sfm_data.IsPoseAndIntrinsicDefined(view))
								continue;
							const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
							const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
							const Vec2 pt = obs_it.second.x;
							bearing.emplace_back((*cam)(cam->get_ud_pixel(pt)));
							poses.emplace_back(pose.asMatrix());
						}
						const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
						Vec4 Xhomogeneous;
						TriangulateNView(bearing_matrix, poses, &Xhomogeneous);
						const Vec3 X = Xhomogeneous.hnormalized();
						// Test validity of the hypothesis (front of the cameras):
						bool bChierality = true;
						int i(0);
						double reprojection_error_sum(0.0);
						for (const auto & obs_it : obs)
						{
							const View * view = sfm_data.views.at(obs_it.first).get();
							if (!sfm_data.IsPoseAndIntrinsicDefined(view))
								continue;

							const Pose3 pose = sfm_data.GetPoseOrDie(view);
							bChierality &= CheiralityTest(bearing[i], pose, X);
							const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
							const Vec2 pt = obs_it.second.x;
							const Vec2 residual = cam->residual(pose(X), pt);
							reprojection_error_sum += residual.norm();
							++i;
						}
						if (bChierality) // Keep the point only if it has a positive depth
						{
							vec_triangulated[landmarkid_it] = X;
							vec_control_points[landmarkid_it] = landmark.X;
							vec_triangulation_errors[landmarkid_it] = reprojection_error_sum / (double)bearing.size();
						}
						else
						{
							std::cout << "Control Point cannot be triangulated (not in front of the cameras)" << std::endl;
							continue;
						}
					}
					if (vec_control_points.size() < 3)
					{
						std::cout << "Insufficient number of triangulated control points." << std::endl;
						std::cout << "Coordinate_id=" << coors_it.first << ",not enough control point(<3),deserted!";
						continue;
					}

					// data conversion to appropriate container
					Mat x1(3, vec_control_points.size()),
						x2(3, vec_control_points.size());
					int col_i(0);
					for (auto & cp_it : vec_control_points)
					{
						x1.col(col_i) = vec_triangulated[cp_it.first];
						x2.col(col_i) = vec_control_points[cp_it.first];
						++col_i;
					}

					std::cout
						<< "Control points observation triangulations:\n"
						<< x1 << std::endl << std::endl
						<< "Control points coords:\n"
						<< x2 << std::endl << std::endl;
					Vec3 t;
					Mat3 R;
					double S;
					if (openMVG::geometry::FindRTS(x2, x1, &S, &t, &R))
					{
						openMVG::geometry::Refine_RTS(x2, x1, &S, &t, &R);
						std::cout << "Found transform:\n"
							<< " scale: " << S << "\n"
							<< " rotation:\n" << R << "\n"
							<< " translation: " << t.transpose() << std::endl;

						const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t / S), S);
						
						seven_parameters[coors_it.first] = sim;

						//如果不把七参数加入BA，我们可以在此直接把控制点坐标换掉
						/*
						std::stringstream os;
						for (auto & cps_it : kzd_coords[coors_it.first])
						{
							const IndexT CPIndex = cps_it;
							// If the control point has not been used, continue...
							if (vec_triangulation_errors.find(CPIndex) == vec_triangulation_errors.end())
								continue;
							Landmark& cp = sfm_data.control_points[cps_it];
							cp.X = sim(vec_control_points[cps_it]);
							os
								<< "CP index: " << CPIndex << "\n"
								<< "CP triangulation error: " << vec_triangulation_errors[CPIndex] << " pixel(s)\n"
								<< "CP registration error: "
								<< (cp.X-vec_triangulated[CPIndex]).norm() << " user unit(s)" << "\n\n";
						}
						std::cout << os.str();	
						*/
					}
				}
			}
		}

	}

	//单坐标系控制点
	if (options.control_point_opt.bUse_control_points&&!options.control_point_opt.bUse_multi_coordinate)
	{
		// Use Ground Control Point:
		//坐标转换
		{
			std::map<IndexT, Vec3> vec_control_points, vec_triangulated;
			std::map<IndexT, double> vec_triangulation_errors;
			for (const auto & control_point_it : sfm_data.control_points)
			{
				/*if(cp_ids.count(control_point_it.first)==0)
				continue;*/
				const Landmark & landmark = control_point_it.second;
				//Triangulate the observations:
				const Observations & obs = landmark.obs;
				std::vector<Vec3> bearing;
				std::vector<Mat34> poses;
				bearing.reserve(obs.size());
				poses.reserve(obs.size());
				for (const auto & obs_it : obs)
				{
					const View * view = sfm_data.views.at(obs_it.first).get();
					if (!sfm_data.IsPoseAndIntrinsicDefined(view))
						continue;
					const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
					const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
					const Vec2 pt = obs_it.second.x;
					bearing.emplace_back((*cam)(cam->get_ud_pixel(pt)));
					poses.emplace_back(pose.asMatrix());
				}
				const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
				Vec4 Xhomogeneous;
				TriangulateNView(bearing_matrix, poses, &Xhomogeneous);
				const Vec3 X = Xhomogeneous.hnormalized();
				// Test validity of the hypothesis (front of the cameras):
				bool bChierality = true;
				int i(0);
				double reprojection_error_sum(0.0);
				for (const auto & obs_it : obs)
				{
					const View * view = sfm_data.views.at(obs_it.first).get();
					if (!sfm_data.IsPoseAndIntrinsicDefined(view))
						continue;

					const Pose3 pose = sfm_data.GetPoseOrDie(view);
					bChierality &= CheiralityTest(bearing[i], pose, X);
					const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
					const Vec2 pt = obs_it.second.x;
					const Vec2 residual = cam->residual(pose(X), pt);
					reprojection_error_sum += residual.norm();
					++i;
				}
				if (bChierality) // Keep the point only if it has a positive depth
				{
					vec_triangulated[control_point_it.first] = X;
					vec_control_points[control_point_it.first] = landmark.X;
					vec_triangulation_errors[control_point_it.first] = reprojection_error_sum / (double)bearing.size();
				}
				else
				{
					std::cout << "Control Point cannot be triangulated (not in front of the cameras)" << std::endl;
					return false;
				}
			}

			if (vec_control_points.size() < 3)
			{
				std::cout << "Insufficient number of triangulated control points." << std::endl;
				return false;
			}

			// compute the similarity
			{
				// data conversion to appropriate container
				Mat x1(3, vec_control_points.size()),
					x2(3, vec_control_points.size());
				for (size_t i = 0; i < vec_control_points.size(); ++i)
				{
					x1.col(i) = vec_triangulated[i];
					x2.col(i) = vec_control_points[i];
				}

				std::cout
					<< "Control points observation triangulations:\n"
					<< x1 << std::endl << std::endl
					<< "Control points coords:\n"
					<< x2 << std::endl << std::endl;

				Vec3 t;
				Mat3 R;
				double S;
				if (openMVG::geometry::FindRTS(x1, x2, &S, &t, &R))
				{
					openMVG::geometry::Refine_RTS(x1, x2, &S, &t, &R);
					std::cout << "Found transform:\n"
						<< " scale: " << S << "\n"
						<< " rotation:\n" << R << "\n"
						<< " translation: " << t.transpose() << std::endl;
					//--
					// Apply the found transformation as a 3D Similarity transformation matrix // S * R * X + t
					//--

					const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t / S), S);
					openMVG::sfm::ApplySimilarity(sim, sfm_data);
				}
				else
				{
					std::cout << "Registration failed. Please check your Control Points coordinates." << std::endl;
					return false;
				}
				vec_control_points.clear();
				vec_triangulated.clear();
				vec_triangulation_errors.clear();
			}
		}
	}

	ceres::Problem problem;

	// Data wrapper for refinement:
	Hash_Map<IndexT, std::vector<double>> map_intrinsics;
	Hash_Map<IndexT, std::vector<double>> map_calibration;
	Hash_Map<IndexT, std::vector<double>> map_sevenparameters;
	Hash_Map<IndexT, std::vector<double>> map_poses;

	//Setup Calibrations &subparametrization-->R_00和T_00固定值，不优化
	//相对位置按照00,01,02,03,04,05,06进行存储
	//0表示主相机，具体对应哪一组view由option中的main_cam决定
	std::cout << "Begin to do the BA adjustment..........." << std::endl;
	std::cout << "The main cam is " << options.main_cam << std::endl;
	std::cout << "The initial values of calibrations:" << std::endl;
	for (int i = 0; i < calibrations.size(); i++)
	{
		const IndexT indexCali = i;
		const Mat3 cali_R = calibrations[i].rotation();
		const Vec3 cali_T = calibrations[i].center();

		double cali_angleAxis[3];
		ceres::RotationMatrixToAngleAxis((const double*)cali_R.data(), cali_angleAxis);
		map_calibration[indexCali] = { cali_angleAxis[0],cali_angleAxis[1],cali_angleAxis[2],cali_T(0),cali_T(1),cali_T(2) };
		std::cout << "Cam" << indexCali << ":" << cali_angleAxis[0] << " " << cali_angleAxis[1] << " " << cali_angleAxis[2] << " " << cali_T(0) << " " << cali_T(1) << " " << cali_T(2) << std::endl;
		double *parameter_block = &map_calibration.at(indexCali)[0];
		problem.AddParameterBlock(parameter_block, 6);
	}
	//R_00,t_00固定
	problem.SetParameterBlockConstant(&map_calibration.at(options.main_cam)[0]);

	//如果将七参数加入BA
	
	std::cout << "Setup Seven parameters data & subparameterization..." << std::endl;
	for (auto & seven_it : seven_parameters)
	{
		const IndexT indexSeven = seven_it.first;
		const Mat3 seven_R = seven_it.second.pose_.rotation();
		const Vec3 seven_T = seven_it.second.pose_.center();

		double seven_angleAxis[3];
		ceres::RotationMatrixToAngleAxis((const double*)seven_R.data(), seven_angleAxis);
		map_sevenparameters[indexSeven] = { seven_angleAxis[0],seven_angleAxis[1],seven_angleAxis[2],seven_T(0),seven_T(1),seven_T(2) };
		std::cout << "Cam" << indexSeven << ":" << seven_angleAxis[0] << " " << seven_angleAxis[1] << " " << seven_angleAxis[2] << " " << seven_T(0) << " " << seven_T(1) << " " << seven_T(2) << std::endl;
		double *parameter_block = &map_sevenparameters.at(indexSeven)[0];
		problem.AddParameterBlock(parameter_block, 6);
		
	}
	//参考坐标系的参数设为常数
	std::cout << "Reference coordinate is :" << reference_coors << std::endl;
	problem.SetParameterBlockConstant(&map_sevenparameters.at(reference_coors)[0]);
	

	std::cout << "Setup Poses data & subparametrization......" << std::endl;
	for (const auto & pose_it : sfm_data.poses)
	{
		//只将相机0的pose参与BA参数优化
		const IndexT indexPose = pose_it.first;
		if (options.IsGrouped)
		{
			if ((int)indexPose%calibrations.size() != options.main_cam)
				continue;
		}
		else
		{
			if ((int)indexPose / station_num != options.main_cam)
				continue;
		}

		const Pose3 & pose = pose_it.second;
		const Mat3 R = pose.rotation();
		const Vec3 t = pose.translation();

		double angleAxis[3];
		ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
		// angleAxis + translation
		std::cout << "indexPose=" << indexPose << ",";
		map_poses[indexPose] = { angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2) };

		double * parameter_block = &map_poses.at(indexPose)[0];
		problem.AddParameterBlock(parameter_block, 6);
		if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
		{
			// set the whole parameter block as constant for best performance
			problem.SetParameterBlockConstant(parameter_block);
		}
		else  // Subset parametrization
		{
			std::vector<int> vec_constant_extrinsic;
			// If we adjust only the translation, we must set ROTATION as constant
			if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
			{
				// Subset rotation parametrization
				vec_constant_extrinsic.insert(vec_constant_extrinsic.end(), { 0,1,2 });
			}
			// If we adjust only the rotation, we must set TRANSLATION as constant
			if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
			{
				// Subset translation parametrization
				vec_constant_extrinsic.insert(vec_constant_extrinsic.end(), { 3,4,5 });
			}
			if (!vec_constant_extrinsic.empty())
			{
				ceres::SubsetParameterization *subset_parameterization =
					new ceres::SubsetParameterization(6, vec_constant_extrinsic);
				problem.SetParameterization(parameter_block, subset_parameterization);
			}
		}
	}

	std::cout<<std::endl;
	std::cout << " Setup Intrinsics data & subparametrization..." << std::endl;

	for (const auto & intrinsic_it : sfm_data.intrinsics)
	{
		const IndexT indexCam = intrinsic_it.first;

		if (isValid(intrinsic_it.second->getType()))
		{
			map_intrinsics[indexCam] = intrinsic_it.second->getParams();
			if (!map_intrinsics.at(indexCam).empty())
			{
				double * parameter_block = &map_intrinsics.at(indexCam)[0];
				problem.AddParameterBlock(parameter_block, map_intrinsics.at(indexCam).size());
				if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
				{
					// set the whole parameter block as constant for best performance
					problem.SetParameterBlockConstant(parameter_block);
				}
				else
				{
					const std::vector<int> vec_constant_intrinsic =
						intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
					if (!vec_constant_intrinsic.empty())
					{
						ceres::SubsetParameterization *subset_parameterization =
							new ceres::SubsetParameterization(
								map_intrinsics.at(indexCam).size(), vec_constant_intrinsic);
						problem.SetParameterization(parameter_block, subset_parameterization);
					}
				}
			}
		}
		else
		{
			std::cerr << "Unsupported camera type." << std::endl;
		}
	}

	// Set a LossFunction to be less penalized by false measurements
	//  - set it to nullptr if you don't want use a lossFunction.
	ceres::LossFunction * p_LossFunction =
		ceres_options_.bUse_loss_function_ ?
		new ceres::HuberLoss(Square(4.0))
		: nullptr;

	if (options.control_point_opt.bUse_control_points)
	{
		std::cout << " Ground Control Point:";
		// - fixed 3D points with weighted observations
		for (auto & gcp_landmark_it : sfm_data.control_points)
		{
			std::set<int> cp_indvec = { 2,9,16,23,30,39,44 };
			if (cp_indvec.count(gcp_landmark_it.first))
				continue;
			std::cout << gcp_landmark_it.first << std::endl;
			const Observations & obs = gcp_landmark_it.second.obs;

			for (const auto & obs_it : obs)
			{
				// Build the residual block corresponding to the track observation:
				const View * view = sfm_data.views.at(obs_it.first).get();
				int cam_index, maincam_pose_index;
				if (options.IsGrouped)
				{
					cam_index = (int)view->id_view % calibrations.size();
					maincam_pose_index = (int)view->id_view / calibrations.size() + options.main_cam;
				}
				else
				{
					cam_index = (int)view->id_view / station_num;
					maincam_pose_index = (int)view->id_view %station_num + station_num*options.main_cam;
				}

				// Each Residual block takes a point and a camera as input and outputs a 2
				// dimensional residual. Internally, the cost function stores the observed
				// image location and compares the reprojection against the observation.
				ceres::CostFunction* cost_function =
					MaincamIntrinsicsToCostFunction(
						sfm_data.intrinsics.at(view->id_intrinsic).get(),
						obs_it.second.x,
						options.control_point_opt.weight);

				if (cost_function)
				{
					//控制点对应的七参与其观测值有关
					if (!map_intrinsics.at(view->id_intrinsic).empty())
					{
						problem.AddResidualBlock(cost_function,
							p_LossFunction,
							&map_intrinsics.at(view->id_intrinsic)[0],
							&map_calibration.at(cam_index)[0],
							&map_sevenparameters.at(obs_it.second.id_feat)[0],
							&map_poses.at(maincam_pose_index)[0],
							gcp_landmark_it.second.X.data());
					}
					else
					{
						std::cerr << "Cannot create a CostFunction for this camera model." << std::endl;
						return false;
					}
				}
			}
			if (obs.empty())
			{
				std::cerr
					<< "Cannot use this GCP id: " << gcp_landmark_it.first
					<< ". There is not linked image observation." << std::endl;
			}
			else
			{
				// Set the 3D point as FIXED (it's a valid GCP)
				problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
			}
		}
	}

	// For all visibility add reprojections errors:
	for (auto & structure_landmark_it : sfm_data.structure)
	{
		const Observations & obs = structure_landmark_it.second.obs;

		for (const auto & obs_it : obs)
		{
			// Build the residual block corresponding to the track observation:
			const View * view = sfm_data.views.at(obs_it.first).get();
			int cam_index, maincam_pose_index;
			if (options.IsGrouped)
			{
				cam_index = (int)view->id_view % calibrations.size();
				maincam_pose_index = (int)view->id_view / calibrations.size() + options.main_cam;
			}
			else
			{
				cam_index = (int)view->id_view / station_num;
				maincam_pose_index = (int)view->id_view %station_num + station_num*options.main_cam;
			}

			// Each Residual block takes a point and a camera as input and outputs a 2
			// dimensional residual. Internally, the cost function stores the observed
			// image location and compares the reprojection against the observation.
			ceres::CostFunction* cost_function =
				MaincamIntrinsicsToCostFunction(sfm_data.intrinsics.at(view->id_intrinsic).get(),
					obs_it.second.x);

			//所有结构点对应的七参数为常数的参考坐标系
			if (cost_function)
			{
				problem.AddResidualBlock(cost_function,
					p_LossFunction,
					&map_intrinsics.at(view->id_intrinsic)[0],
					&map_calibration.at(cam_index)[0],
					&map_sevenparameters.at(reference_coors)[0],
					&map_poses.at(maincam_pose_index)[0],
					structure_landmark_it.second.X.data());
			}
			else
			{
				std::cerr << "Cannot create a CostFunction for this camera model." << std::endl;
				return false;
			}
		}
		if (options.structure_opt == Structure_Parameter_Type::NONE)
			problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());

	}
	

	/*GPS先验
	// Add Pose prior constraints if any
	if (b_usable_prior)
	{
		for (const auto & view_it : sfm_data.GetViews())
		{
			const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
			if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
			{
				// Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
				Vec3 weight;
				//认为XYZ的权值一样
				for (const auto & it_w : Weight_XYZ)
					if (it_w.first == prior->id_pose)
						weight = it_w.second[0] * prior->center_weight_;
				ceres::CostFunction * cost_function =
					new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
						new PoseCenterConstraintCostFunction(prior->pose_center_, weight));

				problem.AddResidualBlock(
					cost_function,
					new ceres::HuberLoss(
						Square(pose_center_robust_fitting_error)),
					&map_poses.at(prior->id_view)[0]);
			}
		}
	}
	*/
	

	// Configure a BA engine and run it
	//  Make Ceres automatically detect the bundle structure.
	ceres::Solver::Options ceres_config_options;
	ceres_config_options.max_num_iterations = 1000;
	ceres_config_options.preconditioner_type =
		static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
	ceres_config_options.linear_solver_type =
		static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
	ceres_config_options.sparse_linear_algebra_library_type =
		static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
	ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
	ceres_config_options.logging_type = ceres::SILENT;
	ceres_config_options.num_threads = ceres_options_.nb_threads_;
#if CERES_VERSION_MAJOR < 2
	ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
#endif
	ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

	// Solve BA
	ceres::Solver::Summary summary;
	ceres::Solve(ceres_config_options, &problem, &summary);
	if (ceres_options_.bCeres_summary_)
		std::cout << summary.FullReport() << std::endl;

	// If no error, get back refined parameters
	if (!summary.IsSolutionUsable())
	{
		if (ceres_options_.bVerbose_)
			std::cout << "Bundle Adjustment failed." << std::endl;
		return false;
	}
	else // Solution is usable
	{
		if (ceres_options_.bVerbose_)
		{
			// Display statistics about the minimization
			std::cout << std::endl
				<< "Bundle Adjustment statistics (approximated RMSE):\n"
				<< " #views: " << sfm_data.views.size() << "\n"
				<< " #poses: " << sfm_data.poses.size() << "\n"
				<< " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
				<< " #tracks: " << sfm_data.structure.size() << "\n"
				<< " #residuals: " << summary.num_residuals << "\n"
				<< " Initial RMSE: " << std::sqrt(summary.initial_cost / summary.num_residuals) << "\n"
				<< " Final RMSE: " << std::sqrt(summary.final_cost / summary.num_residuals) << "\n"
				<< " Time (s): " << summary.total_time_in_seconds << "\n"
				<< std::endl;

			//Display the result of translation calibration refinement
			for (int i=0;i<calibrations.size();i++)
			{
				Vec3 t_refined_cali(map_calibration.at(i)[3], map_calibration.at(i)[4], map_calibration.at(i)[5]);
				std::cout << "The variance of translation:" << t_refined_cali - calibrations[i].translation() << std::endl;
			}
			std::cout << "--------------Display the result of seven parameters----------------" << std::endl;
			for (auto & result_it : map_sevenparameters)
			{
				double * sangleaxis= &result_it.second[0];
				Eigen::Matrix<double, 3, 1> strans;
				strans<<result_it.second[3], result_it.second[4], result_it.second[5];
				Eigen::Matrix<double, 3, 3> srot;
				ceres::AngleAxisToRotationMatrix(sangleaxis, srot.data());
				std::cout << "Rotation:" << "\n"
					<< srot(0, 0) << "," << srot(0, 1) <<","<< srot(0, 2) << "\n"
					<< srot(1, 0) << "," << srot(1, 1) <<","<< srot(1, 2) << "\n"
					<< srot(2, 0) << "," << srot(2, 1) <<","<< srot(2, 2) << std::endl;
				std::cout << "Translation:"
					<< -(srot*strans)(0,0) << "," << -(srot*strans)(1, 0) << "," << -(srot*strans)(2, 0) << std::endl;
				seven_parameters[result_it.first] = Similarity3(Pose3(srot, strans), 1);
			}

			if (options.use_motion_priors_opt)
				std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
		}

		// Update camera poses with refined data
		if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
		{
			for (auto & pose_it : sfm_data.poses)
			{
				//只将相机0的pose参与BA参数优化
				const IndexT indexPose = pose_it.first;
				if (map_poses.count(indexPose))
				{
					Mat3 R_refined;
					ceres::AngleAxisToRotationMatrix(&map_poses.at(indexPose)[0], R_refined.data());
					Vec3 t_refined(map_poses.at(indexPose)[3], map_poses.at(indexPose)[4], map_poses.at(indexPose)[5]);
					// Update the pose
					Pose3 & pose = pose_it.second;
					pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
				}
				else
				{
					//R_cami_j=R_cali_i*R_maincam_j,T_cami_j=T_cali_i+T_miancam_j
					int cam_index, maincam_pose_index;
					if (options.IsGrouped)
					{
						cam_index = (int)indexPose % calibrations.size();
						maincam_pose_index = (int)indexPose / calibrations.size() + options.main_cam;
					}
					else
					{
						cam_index = (int)indexPose / station_num;
						maincam_pose_index = (int)indexPose %station_num + station_num*options.main_cam;
					}
					Mat3 R_refined_maincam;
					ceres::AngleAxisToRotationMatrix(&map_poses.at(maincam_pose_index)[0], R_refined_maincam.data());
					Mat3 R_refined_cali;
					ceres::AngleAxisToRotationMatrix(&map_calibration.at(cam_index)[0], R_refined_cali.data());
					Vec3 t_refined_maincam(map_poses.at(maincam_pose_index)[3], map_poses.at(maincam_pose_index)[4], map_poses.at(maincam_pose_index)[5]);
					Vec3 t_refined_cali(map_calibration.at(cam_index)[3], map_calibration.at(cam_index)[4], map_calibration.at(cam_index)[5]);
					// Update the pose
					Pose3 & pose = pose_it.second;
					pose = Pose3(R_refined_cali*R_refined_maincam, -R_refined_maincam.transpose() * t_refined_maincam + t_refined_cali);
				}

			}
		}

		// Update camera intrinsics with refined data
		if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
		{
			for (auto & intrinsic_it : sfm_data.intrinsics)
			{
				const IndexT indexCam = intrinsic_it.first;

				const std::vector<double> & vec_params = map_intrinsics.at(indexCam);
				intrinsic_it.second->updateFromParams(vec_params);
			}
		}

		// Structure is already updated directly if needed (no data wrapping)
		
		//检查精度
		std::vector<double> vec_residualErrorsX;
		std::vector<double> vec_residualErrorsY;
		std::vector<double> vec_residualErrorsZ;
		for (auto & cp_it : sfm_data.control_points)
		{
			int index_cpt = cp_it.first;
			std::set<int> cp_vec = { 2,9,16,23,30,39,44 };
			if (!cp_vec.count(index_cpt))
				continue;
			openMVG::sfm::Landmark& check_point = sfm_data.control_points[index_cpt];
			std::cout << "the coordinate of the check point(" << index_cpt << "):" << check_point.X[0] << " " << check_point.X[1] << " " << check_point.X[2] << std::endl;
			//Triangulate the observations:
			const openMVG::sfm::Observations & cobs = check_point.obs;
			std::cout << "The number of check_point's visibility:" << cobs.size() << std::endl;
			std::vector<openMVG::Vec3> cbearing;
			std::vector<openMVG::Mat34> cposes;
			double reprojection_error_sum(0);
			cbearing.reserve(cobs.size());
			cposes.reserve(cobs.size());
			for (const auto & cobs_it : cobs)
			{
				const openMVG::sfm::View * view = sfm_data.views.at(cobs_it.first).get();
				std::cout << "View ID:" << view->id_view << "  ";
				if (!sfm_data.IsPoseAndIntrinsicDefined(view))
					continue;
				const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
				const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
				const openMVG::Vec2 pt = cobs_it.second.x;
				std::cout << pt[0] << " " << pt[1];
				cbearing.emplace_back((*cam)(cam->get_ud_pixel(pt)));
				cposes.emplace_back(pose.asMatrix());
				std::cout << std::endl;
			}
			const Eigen::Map<const openMVG::Mat3X> cbearing_matrix(cbearing[0].data(), 3, cbearing.size());
			openMVG::Vec4 cXhomogeneous;

			openMVG::TriangulateNView(cbearing_matrix, cposes, &cXhomogeneous);
			
			const openMVG::Vec3 cX = cXhomogeneous.hnormalized();

			//三角化的重投影误差
			int index_coord = 0;	
			for (const auto & cobs_it : cobs)
			{
				index_coord = cobs_it.second.id_feat;
				const openMVG::sfm::View * view = sfm_data.views.at(cobs_it.first).get();
				if (!sfm_data.IsPoseAndIntrinsicDefined(view))
					continue;
				const openMVG::cameras::IntrinsicBase * cam_ = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
				const openMVG::geometry::Pose3 pose_ = sfm_data.GetPoseOrDie(view);
				const openMVG::Vec2 pt_ = cobs_it.second.x;
				const Vec2 residual = cam_->residual(pose_(cX), pt_);
				reprojection_error_sum += residual.norm();
			}
			double mean_reprojection_error = reprojection_error_sum / (double)cobs.size();
			//对于控制点要进行坐标系转换

			Similarity3 sim = seven_parameters[index_coord];
			openMVG::Vec3 dX = cX - sim(check_point.X);
			std::cout << "the triangulated coordinate of the check point:" << cX[0] << " " << cX[1] << " " << cX[2] << std::endl;
			std::cout << "the mean reprojection error after BA:" << mean_reprojection_error << std::endl;
			std::cout << "the bias after BA:" << dX[0] << " " << dX[1] << " " << dX[2] << std::endl;
			vec_residualErrorsX.push_back(abs(dX[0]));
			vec_residualErrorsY.push_back(abs(dX[1]));
			vec_residualErrorsZ.push_back(abs(dX[2]));
		}
		std::cout << std::endl << "\nReprojectX error statistics : \n ";
		minMaxMeanMedian<double>(vec_residualErrorsX.begin(), vec_residualErrorsX.end());

		std::cout << std::endl << "\nReprojectY error statistics : \n ";
		minMaxMeanMedian<double>(vec_residualErrorsY.begin(), vec_residualErrorsY.end());

		std::cout << std::endl << "\nReprojectZ error statistics : \n ";
		minMaxMeanMedian<double>(vec_residualErrorsZ.begin(), vec_residualErrorsZ.end());
		if (b_usable_prior)
		{
			// set back to the original scene centroid
			openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

			//--
			// - Compute some fitting statistics
			//--

			// Collect corresponding camera centers
			std::vector<Vec3> X_SfM, X_GPS;
			for (const auto & view_it : sfm_data.GetViews())
			{
				const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
				if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
				{
					X_SfM.push_back(sfm_data.GetPoses().at(prior->id_pose).center());
					X_GPS.push_back(prior->pose_center_);
				}
			}
			// Compute the registration fitting error (once BA with Prior have been used):
			if (X_GPS.size() > 3)
			{
				// Compute the median residual error
				Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
				std::cout
					<< "Pose prior statistics (user units):\n"
					<< " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
					<< " - Final fitting error:";
				minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
			}
		}
		return true;
	}
}

void Bundle_Adjustment_Ceres::Check_residual(sfm::SfM_Data & sfm_data,std::vector<geometry::Pose3>& calibrations)
{
	//构建所有参数
	Hash_Map<IndexT, std::vector<double>> map_intrinsics;
	Hash_Map<IndexT, std::vector<double>> map_calibration;
	Hash_Map<IndexT, std::vector<double>> map_sevens;
	Hash_Map<IndexT, std::vector<double>> map_poses;
	
	for (int i = 0; i < calibrations.size(); i++)
	{
		const IndexT indexCali = i;
		const Mat3 cali_R = calibrations[i].rotation();
		const Vec3 cali_T = calibrations[i].center();

		double cali_angleAxis[3];
		ceres::RotationMatrixToAngleAxis((const double*)cali_R.data(), cali_angleAxis);
		map_calibration[indexCali] = { cali_angleAxis[0],cali_angleAxis[1],cali_angleAxis[2],cali_T(0),cali_T(1),cali_T(2) };
		std::cout << "Cam" << indexCali << ":" << cali_angleAxis[0] << " " << cali_angleAxis[1] << " " << cali_angleAxis[2] << " " << cali_T(0) << " " << cali_T(1) << " " << cali_T(2) << std::endl;
	}
	
	//七参数估计
	{
		//控制点为多坐标系
		std::map<int, std::set<int>>kzd_coords;
		std::map<int, geometry::Similarity3> seven_parameters;
		for (auto & gcp_landmark_it : sfm_data.control_points)
		{
			const Observations & obs = gcp_landmark_it.second.obs;
			kzd_coords[obs.begin()->second.id_feat].insert(gcp_landmark_it.first);
		}

		//统计每个坐标系的观测值数
		std::map<int, int>coords_obs;
		for (auto coors_it = kzd_coords.begin(); coors_it != kzd_coords.end(); ++coors_it)
		{
			//删除控制点数目少于3的坐标系
			if (coors_it->second.size() < 3)
			{
				kzd_coords.erase(coors_it);
				continue;
			}
			int count = 0;
			//删除观测值少于2的控制点
			for (auto landmarkid_it : coors_it->second)
			{
				int obs_size = sfm_data.control_points.at(landmarkid_it).obs.size();
				if (obs_size < 2)
				{
					coors_it->second.erase(landmarkid_it);
					continue;
				}
				count += obs_size;
			}

			coords_obs[coors_it->first] = count;
			std::cout << "The coords_" << coors_it->first << " ,obs:" << coords_obs[coors_it->first] << ",control_points_num:" << coors_it->second.size() << std::endl;
		}

		//然后以控制点平均观测值最多的一个坐标系为参考
		int reference_coors = 0;
		double avg_obs_max = coords_obs[reference_coors] / kzd_coords[reference_coors].size();
		for (auto & refer_it : coords_obs)
		{
			double avg_obs = 1.0*refer_it.second / kzd_coords[refer_it.first].size();
			std::cout << "The coords_" << refer_it.first << " has average obs:" << avg_obs << std::endl;
			if (avg_obs > avg_obs_max)
				reference_coors = refer_it.first;
		}
		std::cout << "The reference coordinate is:" << reference_coors << std::endl;

		//参考坐标系的七参数
		seven_parameters[reference_coors] = geometry::Similarity3();

		//进行七参转换
		{
			std::cout << "..............................The first transforming..............................." << std::endl;
			std::map<IndexT, Vec3> vec_control_points, vec_triangulated;
			std::map<IndexT, double> vec_triangulation_errors;
			for (auto & landmarkid_it : kzd_coords[reference_coors])
			{
				std::cout << "landmarkid_it=" << landmarkid_it << std::endl;
				const openMVG::sfm::Landmark& landmark = sfm_data.control_points[landmarkid_it];
				const Observations & obs = landmark.obs;
				std::cout << "Triangulate the observations(" << obs.size() << ").................." << std::endl;
				std::vector<Vec3> bearing;
				std::vector<Mat34> poses;
				bearing.reserve(obs.size());
				poses.reserve(obs.size());
				for (const auto & obs_it : obs)
				{
					const View * view = sfm_data.views.at(obs_it.first).get();
					if (!sfm_data.IsPoseAndIntrinsicDefined(view))
						continue;
					const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
					const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
					const Vec2 pt = obs_it.second.x;
					bearing.emplace_back((*cam)(cam->get_ud_pixel(pt)));
					poses.emplace_back(pose.asMatrix());
				}
				const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
				Vec4 Xhomogeneous;
				TriangulateNView(bearing_matrix, poses, &Xhomogeneous);
				const Vec3 X = Xhomogeneous.hnormalized();
				// Test validity of the hypothesis (front of the cameras):
				bool bChierality = true;
				int i(0);
				double reprojection_error_sum(0.0);
				for (const auto & obs_it : obs)
				{
					const View * view = sfm_data.views.at(obs_it.first).get();
					if (!sfm_data.IsPoseAndIntrinsicDefined(view))
						continue;

					const Pose3 pose = sfm_data.GetPoseOrDie(view);
					bChierality &= CheiralityTest(bearing[i], pose, X);
					const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
					const Vec2 pt = obs_it.second.x;
					const Vec2 residual = cam->residual(pose(X), pt);
					reprojection_error_sum += residual.norm();
					++i;
				}
				if (bChierality) // Keep the point only if it has a positive depth
				{
					vec_triangulated[landmarkid_it] = X;
					std::cout << "vec_triangulated[" << landmarkid_it << "] =" << X[0] << "," << X[1] << "," << X[2] << std::endl;
					vec_control_points[landmarkid_it] = landmark.X;
					std::cout << "vec_control_points[" << landmarkid_it << "] =" << landmark.X[0] << "," << landmark.X[1] << "," << landmark.X[2] << std::endl;
					vec_triangulation_errors[landmarkid_it] = reprojection_error_sum / (double)bearing.size();
				}
				else
				{
					std::cout << "Control Point cannot be triangulated (not in front of the cameras)" << std::endl;
					return;
				}
			}
			if (vec_control_points.size() < 3)
			{
				std::cout << "Insufficient number of triangulated control points." << std::endl;
				return;
			}

			// data conversion to appropriate container
			Mat x1(3, vec_control_points.size()),
				x2(3, vec_control_points.size());
			int col_i(0);
			for (auto & cp_it : vec_control_points)
			{
				x1.col(col_i) = vec_triangulated[cp_it.first];
				x2.col(col_i) = vec_control_points[cp_it.first];
				++col_i;
			}

			std::cout
				<< "Control points observation triangulations:\n"
				<< x1 << std::endl << std::endl
				<< "Control points coords:\n"
				<< x2 << std::endl << std::endl;
			Vec3 t;
			Mat3 R;
			double S;
			if (openMVG::geometry::FindRTS(x1, x2, &S, &t, &R))
			{
				openMVG::geometry::Refine_RTS(x1, x2, &S, &t, &R);
				std::cout << "Found transform:\n"
					<< " scale: " << S << "\n"
					<< " rotation:\n" << R << "\n"
					<< " translation: " << t.transpose() << std::endl;
				//--
				// Apply the found transformation as a 3D Similarity transformation matrix // S * R * X + t
				//--

				const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t / S), S);
				openMVG::sfm::ApplySimilarity(sim, sfm_data);

				// Display some statistics:
				std::stringstream os;
				for (auto & landmarkid_it : kzd_coords[reference_coors])
				{
					const IndexT CPIndex = landmarkid_it;
					// If the control point has not been used, continue...
					if (vec_triangulation_errors.find(CPIndex) == vec_triangulation_errors.end())
						continue;

					os
						<< "CP index: " << CPIndex << "\n"
						<< "CP triangulation error: " << vec_triangulation_errors[CPIndex] << " pixel(s)\n"
						<< "CP registration error: "
						<< (sim(vec_triangulated[CPIndex]) - vec_control_points[CPIndex]).norm() << " user unit(s)" << "\n\n";
				}
				std::cout << os.str();
			}
		}

		//估计其他七参数
		for (auto & coors_it : kzd_coords)
		{
			std::cout << "!!!!!!!!!!!!!!Beginning to evaluate the seven parameters of Coordinate_" << coors_it.first << std::endl;
			if (coors_it.first == reference_coors)
				continue;
			//控制点数少于三的坐标系不参与计算
			if (coors_it.second.size() < 3)
			{
				std::cout << "Coordinate_id=" << coors_it.first << ",not enough control point(<3),deserted!";
				continue;
			}
			//七参计算
			{
				std::map<IndexT, Vec3> vec_control_points, vec_triangulated;
				std::map<IndexT, double> vec_triangulation_errors;
				for (auto & landmarkid_it : kzd_coords[coors_it.first])
				{
					const openMVG::sfm::Landmark& landmark = sfm_data.control_points[landmarkid_it];
					//Triangulate the observations:
					const Observations & obs = landmark.obs;
					std::vector<Vec3> bearing;
					std::vector<Mat34> poses;
					bearing.reserve(obs.size());
					poses.reserve(obs.size());
					for (const auto & obs_it : obs)
					{
						const View * view = sfm_data.views.at(obs_it.first).get();
						if (!sfm_data.IsPoseAndIntrinsicDefined(view))
							continue;
						const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
						const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
						const Vec2 pt = obs_it.second.x;
						bearing.emplace_back((*cam)(cam->get_ud_pixel(pt)));
						poses.emplace_back(pose.asMatrix());
					}
					const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
					Vec4 Xhomogeneous;
					TriangulateNView(bearing_matrix, poses, &Xhomogeneous);
					const Vec3 X = Xhomogeneous.hnormalized();
					// Test validity of the hypothesis (front of the cameras):
					bool bChierality = true;
					int i(0);
					double reprojection_error_sum(0.0);
					for (const auto & obs_it : obs)
					{
						const View * view = sfm_data.views.at(obs_it.first).get();
						if (!sfm_data.IsPoseAndIntrinsicDefined(view))
							continue;

						const Pose3 pose = sfm_data.GetPoseOrDie(view);
						bChierality &= CheiralityTest(bearing[i], pose, X);
						const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
						const Vec2 pt = obs_it.second.x;
						const Vec2 residual = cam->residual(pose(X), pt);
						reprojection_error_sum += residual.norm();
						++i;
					}
					if (bChierality) // Keep the point only if it has a positive depth
					{
						vec_triangulated[landmarkid_it] = X;
						vec_control_points[landmarkid_it] = landmark.X;
						vec_triangulation_errors[landmarkid_it] = reprojection_error_sum / (double)bearing.size();
					}
					else
					{
						std::cout << "Control Point cannot be triangulated (not in front of the cameras)" << std::endl;
						continue;
					}
				}
				if (vec_control_points.size() < 3)
				{
					std::cout << "Insufficient number of triangulated control points." << std::endl;
					std::cout << "Coordinate_id=" << coors_it.first << ",not enough control point(<3),deserted!";
					continue;
				}

				// data conversion to appropriate container
				Mat x1(3, vec_control_points.size()),
					x2(3, vec_control_points.size());
				int col_i(0);
				for (auto & cp_it : vec_control_points)
				{
					x1.col(col_i) = vec_triangulated[cp_it.first];
					x2.col(col_i) = vec_control_points[cp_it.first];
					++col_i;
				}

				std::cout
					<< "Control points observation triangulations:\n"
					<< x1 << std::endl << std::endl
					<< "Control points coords:\n"
					<< x2 << std::endl << std::endl;
				Vec3 t;
				Mat3 R;
				double S;
				if (openMVG::geometry::FindRTS(x2, x1, &S, &t, &R))
				{
					openMVG::geometry::Refine_RTS(x2, x1, &S, &t, &R);
					std::cout << "Found transform:\n"
						<< " scale: " << S << "\n"
						<< " rotation:\n" << R << "\n"
						<< " translation: " << t.transpose() << std::endl;

					const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t / S), S);

					seven_parameters[coors_it.first] = sim;
				}
			}
		}
		for (auto & seven_it : seven_parameters)
		{
			const IndexT indexSeven = seven_it.first;
			const Mat3 seven_R = seven_it.second.pose_.rotation();
			const Vec3 seven_T = seven_it.second.pose_.center();

			double seven_angleAxis[3];
			ceres::RotationMatrixToAngleAxis((const double*)seven_R.data(), seven_angleAxis);
			map_sevens[indexSeven] = { seven_angleAxis[0],seven_angleAxis[1],seven_angleAxis[2],seven_T(0),seven_T(1),seven_T(2) };
		}
	}

	for (const auto & pose_it : sfm_data.poses)
	{
		const IndexT indexPose = pose_it.first;

		const Pose3 & pose = pose_it.second;
		const Mat3 R = pose.rotation();
		const Vec3 t = pose.translation();

		double angleAxis[3];
		ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
		// angleAxis + translation
		map_poses[indexPose] = { angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2) };
	}

	for (const auto & intrinsic_it : sfm_data.intrinsics)
	{
		const IndexT indexCam = intrinsic_it.first;

		if (isValid(intrinsic_it.second->getType()))
			map_intrinsics[indexCam] = intrinsic_it.second->getParams();
		else
			std::cerr << "Unsupported camera type." << std::endl;
	}

	//计算每个观测值残差
	Landmark & cp = sfm_data.control_points.at(12);
	std::cout << "kzd positionXYZ:" << cp.X[0] << " " << cp.X[1] << " " << cp.X[2] << std::endl;
	for (auto & cpobs_it:cp.obs)
	{
		int cam_id = cpobs_it.first/49;
		int maincam_Index= (int)cpobs_it.first%49 + 49*5;
		const double * cam_cR = &map_calibration.at(cam_id)[0];
		Eigen::Map<const Eigen::Matrix<double, 3, 1>> cam_ct(&map_calibration.at(cam_id)[3]);
		std::cout <<"View id="<<cpobs_it.first<< ",cam_id of this obs:" << cam_id <<std::endl;
		Eigen::Matrix<double, 3, 3> R_i0;
		ceres::AngleAxisToRotationMatrix(cam_cR, R_i0.data());
		const double * cam_R = &map_poses.at(maincam_Index)[0];
		Eigen::Map<const Eigen::Matrix<double, 3, 1>> cam_t(&map_poses.at(maincam_Index)[3]);
		Eigen::Matrix<double, 3, 3> R_0;
		ceres::AngleAxisToRotationMatrix(cam_R, R_0.data());

		Eigen::Matrix<double, 3, 3> R_i = R_i0*R_0;
		double angleAxis[3];
		ceres::RotationMatrixToAngleAxis(R_i.data(), angleAxis);
		std::cout << "angleAxis:" << angleAxis[0] << "," << angleAxis[1] << "," << angleAxis[2] << std::endl;

		double pos_3dpoint[3] = { cp.X[0],cp.X[1],cp.X[2] };
		std::cout << "reference coordinate: " << cpobs_it.second.id_feat << std::endl;
		const double * cam_cS = &map_sevens[cpobs_it.second.id_feat][0];
		const Eigen::Matrix<double, 3, 1> cam_cSt(&map_sevens[cpobs_it.second.id_feat][3]);
		Eigen::Matrix<double, 3, 3> R_S;
		ceres::AngleAxisToRotationMatrix(cam_cS, R_S.data());
		std::cout << "Rotation:" << "\n"
			<< R_S(0, 0) << "," << R_S(0, 1) <<","<< R_S(0, 2) << "\n"
			<< R_S(1, 0) << "," << R_S(1, 1) <<","<< R_S(1, 2) << "\n"
			<< R_S(2, 0) << "," << R_S(2, 1) <<","<< R_S(2, 2) << std::endl;
		std::cout << "Translation:"
			<< -(R_S*cam_cSt)[0] << "," << -(R_S*cam_cSt)[1] << "," << -(R_S*cam_cSt)[2] << std::endl;
		//将pos_3dpoint的坐标系转换
		const Eigen::Matrix<double, 3, 1> pos_3d_point(&pos_3dpoint[0]);
		Eigen::Matrix<double, 3, 1> pos_3dpointS = R_S*(pos_3dpointS - cam_cSt);

		// Rotate the point according the camera rotation
		Eigen::Matrix<double, 3, 1> transformed_point;
		ceres::AngleAxisRotatePoint(angleAxis, pos_3dpointS.data(), transformed_point.data());
		// Apply the camera translation
		transformed_point += R_i*(R_0.transpose()*cam_t - cam_ct);
		std::cout << "transformed_point:" << transformed_point(0, 0) << "," << transformed_point(1,0) << "," << transformed_point(2, 0) << std::endl;
		
		//Calculate the real point from known pose
		Eigen::Matrix<double, 3, 1> transformed_point_known;
		const double * cam_R_known = &map_poses.at(cpobs_it.first)[0];
		std::cout << "cam_R_known:" << cam_R_known[0] << "," << cam_R_known[1] << "," << cam_R_known[2] << std::endl;
		Eigen::Map<const Eigen::Matrix<double, 3, 1>> cam_t_known(&map_poses.at(cpobs_it.first)[3]);
		ceres::AngleAxisRotatePoint(cam_R_known, pos_3dpointS.data(), transformed_point_known.data());
		transformed_point_known += cam_t_known;
		std::cout << "transformed_point_known:" << transformed_point_known(0, 0) << "," << transformed_point_known(1, 0) << "," << transformed_point_known(2, 0) << std::endl;

		// Transform the point from homogeneous to euclidean (undistorted point)
		const Eigen::Matrix<double, 2, 1> projected_point = transformed_point.hnormalized();
		//--
		// Apply intrinsic parameters
		//--
		int index_intrinsic = sfm_data.views.at(cpobs_it.first).get()->id_intrinsic;
		const double& focal = map_intrinsics.at(index_intrinsic)[0];
		const double& principal_point_x = map_intrinsics.at(index_intrinsic)[1];
		const double& principal_point_y = map_intrinsics.at(index_intrinsic)[2];
		const double& k1 = map_intrinsics.at(index_intrinsic)[3];
		const double& k2 = map_intrinsics.at(index_intrinsic)[4];
		const double& k3 = map_intrinsics.at(index_intrinsic)[5];
		const double& k4 = map_intrinsics.at(index_intrinsic)[6];
		//const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
		//const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

		// Apply distortion (xd,yd) = disto(x_u,y_u)
		//鱼眼相机
		const double r2 = projected_point.squaredNorm();
		const double r = sqrt(r2);
		const double
			theta = atan(r),
			theta2 = theta*theta,
			theta3 = theta2*theta,
			theta4 = theta2*theta2,
			theta5 = theta4*theta,
			theta7 = theta3*theta3*theta, //thetha6*theta
			theta8 = theta4*theta4,
			theta9 = theta8*theta;
		const double theta_dist = theta + k1*theta3 + k2*theta5 + k3*theta7 + k4*theta9;
		const double inv_r = r > double(1e-8) ? double(1.0) / r : double(1.0);
		const double cdist = r > double(1e-8) ? theta_dist * inv_r : double(1.0);
		Vec2 residual(principal_point_x + (projected_point.x() * cdist) * focal - cpobs_it.second.x[0],
			principal_point_y + (projected_point.y() * cdist) * focal - cpobs_it.second.x[1]);

		//布朗模型
		/*Vec2 residual( principal_point_x + (projected_point.x() * r_coeff + t_x) * focal - cpobs_it.second.x[0],
		principal_point_y + (projected_point.y() * r_coeff + t_y) * focal - cpobs_it.second.x[1]);*/
		std::cout << "residual:" << residual[0] << "," << residual[1] << std::endl;
		std::cout << std::endl;
	}
	return ;
}

} // namespace sfm
} // namespace openMVG
