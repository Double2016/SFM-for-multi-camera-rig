// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_BA_ceres.hpp"

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

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

#include "Bundle_Adjustment_Ceres_group.h"

#include <ceres/rotation.h>
#include <ceres/types.h>

#include <iostream>
#include <limits>

namespace openMVG {
	namespace sfm {

		using namespace openMVG::cameras;
		using namespace openMVG::geometry;

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
				OFFSET_PRINCIPAL_POINT_Y = 2
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
				const T *cam_cR = cam_calibration;
				Eigen::Map<const Eigen::Matrix<T, 3, 1>> cam_ct(&cam_calibration[3]);
				Mat3 R_i0;
				ceres::AngleAxisToRotationMatrix(cam_cR,Rij.data());

				const T * cam_R = cam_extrinsics;
				Eigen::Map<const Eigen::Matrix<T, 3, 1>> cam_t(&cam_extrinsics[3]);

				Mat3 R_0;
				ceres::AngleAxisToRotationMatrix(cam_R, R_0.data());

				Mat3 R = R_i0*R_0;
				double angleAxis[3];
				ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);

				Eigen::Matrix<T, 3, 1> transformed_point;
				// Rotate the point according the camera rotation
				ceres::AngleAxisRotatePoint(&angleAxis, pos_3dpoint, transformed_point.data());

				// Apply the camera translation
				transformed_point += R*(R_i0*R_0.transpose()*cam_t-*cam_ct);

				// Transform the point from homogeneous to euclidean (undistorted point)
				const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

				//--
				// Apply intrinsic parameters
				//--

				const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
				const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
				const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

				// Apply focal length and principal point to get the final image coordinates

				// Compute and return the error is the difference between the predicted
				//  and observed position
				Eigen::Map<Eigen::Matrix<T, 2, 1>> residuals(out_residuals);
				residuals << principal_point_x + projected_point.x() * focal - m_pos_2dpoint[0],
					principal_point_y + projected_point.y() * focal - m_pos_2dpoint[1];
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
							<ResidualErrorFunctor_Pinhole_Intrinsic, 2, 3, 6, 3>(
								new ResidualErrorFunctor_Pinhole_Intrinsic(observation.data())));
				}
				else
				{
					return
						(new ceres::AutoDiffCostFunction
							<WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic>, 2, 3, 6, 3>
							(new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic>
							(new ResidualErrorFunctor_Pinhole_Intrinsic(observation.data()), weight)));
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
			return ResidualErrorFunctor_Group_Intrinsic::Create(observation, weight);
		}

		Bundle_Adjustment_Ceres_group::BA_Ceres_options::BA_Ceres_options
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


		Bundle_Adjustment_Ceres_group::Bundle_Adjustment_Ceres_group
		(
			const Bundle_Adjustment_Ceres_group::BA_Ceres_options & options
		)
			: ceres_options_(options)
		{}

		Bundle_Adjustment_Ceres_group::BA_Ceres_options &
			Bundle_Adjustment_Ceres_group::ceres_options()
		{
			return ceres_options_;
		}

		bool Bundle_Adjustment_Ceres_group::Adjust
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
			if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
			{
				// - Compute a robust X-Y affine transformation & apply it
				// - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
				{
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
							b_usable_prior = true; // PRIOR can be used safely

												   // Compute the median residual error once the registration is applied
							for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
							{
								pos = sim(pos);
							}
							Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
							std::sort(residual.data(), residual.data() + residual.size());
							pose_center_robust_fitting_error = residual(residual.size() / 2);

							// Apply the found transformation to the SfM Data Scene
							openMVG::sfm::ApplySimilarity(sim, sfm_data);

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

			if (options.control_point_opt.bUse_control_points)
			{
				// Use Ground Control Point:
				// - fixed 3D points with weighted observations
				for (auto & gcp_landmark_it : sfm_data.control_points)
				{
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

			//// Add Pose prior constraints if any
			//if (b_usable_prior)
			//{
			//	for (const auto & view_it : sfm_data.GetViews())
			//	{
			//		const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
			//		if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
			//		{
			//			// Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
			//			ceres::CostFunction * cost_function =
			//				new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
			//					new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

			//			problem.AddResidualBlock(
			//				cost_function,
			//				new ceres::HuberLoss(
			//					Square(pose_center_robust_fitting_error)),
			//				&map_poses.at(prior->id_view)[0]);
			//		}
			//	}
			//}

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
						<< " Initial RMSE: " << std::sqrt(summary.initial_cost / summary.num_residuals) << "\n"
						<< " Final RMSE: " << std::sqrt(summary.final_cost / summary.num_residuals) << "\n"
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

				// Structure is already updated directly if needed (no data wrapping)

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


		bool Bundle_Adjustment_Ceres_group::Adjust_group
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

			ceres::Problem problem;

			// Data wrapper for refinement:
			Hash_Map<IndexT, std::vector<double>> map_intrinsics;
			Hash_Map<IndexT, std::vector<double>> map_calibration;
			Hash_Map<IndexT, std::vector<double>> map_poses;
			
			//Setup Calibrations &subparametrization-->R_00和T_00固定值，不优化
			//相对位置按照00,01,02,03,04,05,06进行存储
			//0表示主相机，具体对应哪一组view由option中的main_cam决定
			for (int i=0;i<calibrations.size();i++)
			{
				const IndexT indexCali = i;
				const Mat3 cali_R = calibrations[i].rotation();
				const Vec3 cali_T = calibrations[i].center();
				double cali_angleAxis[3];
				ceres::RotationMatrixToAngleAxis((const double*)cali_R.data(), cali_angleAxis);
				map_calibration[indexCali] = { cali_angleAxis[0],cali_angleAxis[1],cali_angleAxis[2],cali_T(0),cali_T(1),cali_T(2) };
			
				double *parameter_block = &map_poses.at(indexCali)[0];
				problem.AddParameterBlock(parameter_block, 6);
			}

			// Setup Poses data & subparametrization
			for (const auto & pose_it : sfm_data.poses)
			{
				//只将相机0的pose参与BA参数优化
				const IndexT indexPose = pose_it.first;
				if((int)indexPose%calibrations.size()!=options.main_cam)
					continue;

				const Pose3 & pose = pose_it.second;
				const Mat3 R = pose.rotation();
				const Vec3 t = pose.translation();

				double angleAxis[3];
				ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
				// angleAxis + translation
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

			// Setup Intrinsics data & subparametrization
			//只加入相机main的内参数
			const IndexT indexCam = sfm_data.views.at(options.main_cam).get()->id_intrinsic;
			if (isValid(sfm_data.GetIntrinsics().at(indexCam)->getType()))
			{
				map_intrinsics[indexCam] = sfm_data.GetIntrinsics().at(indexCam)->getParams();
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
							sfm_data.GetIntrinsics().at(indexCam)->subsetParameterization(options.intrinsics_opt);
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

			// Set a LossFunction to be less penalized by false measurements
			//  - set it to nullptr if you don't want use a lossFunction.
			ceres::LossFunction * p_LossFunction =
				ceres_options_.bUse_loss_function_ ?
				new ceres::HuberLoss(Square(4.0))
				: nullptr;

			// For all visibility add reprojections errors:
			for (auto & structure_landmark_it : sfm_data.structure)
			{
				const Observations & obs = structure_landmark_it.second.obs;

				for (const auto & obs_it : obs)
				{
					// Build the residual block corresponding to the track observation:
					const View * view = sfm_data.views.at(obs_it.first).get();
					const int cam_index = (int)view->id_view % calibrations.size();
					int station_index = (int)view->id_view / calibrations.size();

					// Each Residual block takes a point and a camera as input and outputs a 2
					// dimensional residual. Internally, the cost function stores the observed
					// image location and compares the reprojection against the observation.
					ceres::CostFunction* cost_function =
						IntrinsicsToCostFunction(sfm_data.intrinsics.at(view->id_intrinsic).get(),
							obs_it.second.x);

					if (cost_function)
					{
						problem.AddResidualBlock(cost_function,
							p_LossFunction,
							&map_intrinsics.at(view->id_intrinsic)[0],
							&map_calibration.at(cam_index)[0],
							&map_poses.at(station_index*7+options.main_cam)[0],
							structure_landmark_it.second.X.data());
					}
					else
					{
						std::cerr << "Cannot create a CostFunction for this camera model." << std::endl;
						return false;
					}
					//R_00,t_00固定
					if (cam_index==options.main_cam)
						problem.SetParameterBlockConstant(&map_calibration.at(options.main_cam)[0]);
				}
				if (options.structure_opt == Structure_Parameter_Type::NONE)
					problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
				
			}

			if (options.control_point_opt.bUse_control_points)
			{
				// Use Ground Control Point:
				// - fixed 3D points with weighted observations
				for (auto & gcp_landmark_it : sfm_data.control_points)
				{
					const Observations & obs = gcp_landmark_it.second.obs;

					for (const auto & obs_it : obs)
					{
						// Build the residual block corresponding to the track observation:
						const View * view = sfm_data.views.at(obs_it.first).get();
						const int cam_index = (int)view->id_view % calibrations.size();

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
									p_LossFunction,
									&map_intrinsics.at(view->id_intrinsic)[0],
									&map_calibration.at(cam_index)[0],
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

			// Add Pose prior constraints if any
			//if (b_usable_prior)
			//{
			//	for (const auto & view_it : sfm_data.GetViews())
			//	{
			//		const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
			//		if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
			//		{
			//			// Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
			//			Vec3 weight;
			//			//认为XYZ的权值一样
			//			for (const auto & it_w : Weight_XYZ)
			//				if (it_w.first == prior->id_pose)
			//					weight = it_w.second[0] * prior->center_weight_;
			//			ceres::CostFunction * cost_function =
			//				new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
			//					new PoseCenterConstraintCostFunction(prior->pose_center_, weight));

			//			problem.AddResidualBlock(
			//				cost_function,
			//				new ceres::HuberLoss(
			//					Square(pose_center_robust_fitting_error)),
			//				&map_poses.at(prior->id_view)[0]);
			//		}
			//	}
			//}

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
						<< " Initial RMSE: " << std::sqrt(summary.initial_cost / summary.num_residuals) << "\n"
						<< " Final RMSE: " << std::sqrt(summary.final_cost / summary.num_residuals) << "\n"
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

				// Structure is already updated directly if needed (no data wrapping)

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

	} // namespace sfm
} // namespace openMVG
