#pragma once

#include "D:\yss\Empty-three-orientation\src\openMVG\sfm\sfm_data_BA.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/cameras/Camera_Common.hpp"

namespace ceres { class CostFunction; }
namespace openMVG { namespace cameras { struct IntrinsicBase; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {
	namespace sfm {

		/// Create the appropriate cost functor according the provided input camera intrinsic model
		/// Can be residual cost functor can be weighetd if desired (default 0.0 means no weight).
		/*ceres::CostFunction * IntrinsicsToCostFunction
		(
			cameras::IntrinsicBase * intrinsic,
			const Vec2 & observation,
			const double weight = 0.0
		);*/

		class Bundle_Adjustment_Ceres_group :
			public Bundle_Adjustment
		{
		public:
			struct BA_Ceres_options
			{
				bool bVerbose_;
				unsigned int nb_threads_;
				bool bCeres_summary_;
				int linear_solver_type_;
				int preconditioner_type_;
				int sparse_linear_algebra_library_type_;
				double parameter_tolerance_;
				bool bUse_loss_function_;

				BA_Ceres_options(const bool bVerbose = true, bool bmultithreaded = true);
			};
		private:
			BA_Ceres_options ceres_options_;
		public:
			explicit	Bundle_Adjustment_Ceres_group(
				const Bundle_Adjustment_Ceres_group::BA_Ceres_options & options =
				std::move(BA_Ceres_options()));

			BA_Ceres_options & ceres_options();

			bool Adjust_group
			(
				// the SfM scene to refine
				sfm::SfM_Data & sfm_data,
				//the cameras' calibration initial values
				std::vector<geometry::Pose3>& calibrations,
				// tell which parameter needs to be adjusted
				const Optimize_Options & options
			) ;

			bool Adjust
			(
				// the SfM scene to refine
				sfm::SfM_Data & sfm_data,
				// tell which parameter needs to be adjusted
				const Optimize_Options & options
			) override;


			~Bundle_Adjustment_Ceres_group();
		};
	}
}




