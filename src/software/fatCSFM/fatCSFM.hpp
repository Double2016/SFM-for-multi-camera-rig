//sfm_init_ImageListing
#include "openMVG/cameras/cameras.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"
#include "openMVG/exif/sensor_width_database/ParseDatabase.hpp"
#include "openMVG/geodesy/geodesy.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/types.hpp"


#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <utility>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::geodesy;
using namespace openMVG::image;
using namespace openMVG::sfm;

//features
#include <cereal/archives/json.hpp>
#include "openMVG/features/akaze/image_describer_akaze_io.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer_io.hpp"
#include "openMVG/features/regions_factory_io.hpp"
#include "openMVG/system/timer.hpp"
#include "nonFree/sift/SIFT_describer_io.hpp"
#include <cereal/details/helpers.hpp>
#include <atomic>
#include <cstdlib>
#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif
using namespace openMVG::features;
using namespace std;
//q
#include "siftgpu_GLSL.h"

//matches
#include "openMVG/graph/graph.hpp"
#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/matching_image_collection/F_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust_Angular.hpp"
#include "openMVG/matching_image_collection/Eo_Robust.hpp"
#include "openMVG/matching_image_collection/H_ACRobust.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/geometry/Similarity3_Kernel.hpp"
#include "openMVG/robust_estimation/robust_estimator_LMeds.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/geometry/rigid_transformation3D_srt.hpp"

using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace openMVG::matching_image_collection;
//global sfm
#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_report.hpp"

//incrementalSFM
#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"

//incrementalSFM2
#include "openMVG/sfm/pipelines/sequential/sequential_SfM2.hpp"
#include "openMVG/sfm/pipelines/sequential/SfmSceneInitializerMaxPair.hpp"
#include "openMVG/sfm/pipelines/sequential/SfmSceneInitializerStellar.hpp"

//GroupSFM
#include "openMVG/sfm/sfm_data_transform.hpp"

enum EPairMode
{
	PAIR_EXHAUSTIVE = 0,
	PAIR_CONTIGUOUS = 1,
	PAIR_FOR_MULTICAMERAS=2
};


class CSFM
{
public:
    CSFM();
	~CSFM();
	//Global Variants
	SfM_Data sfm_data;
	
	std::string sImageDir,
		sfileDatabase = "",
		sOutputDir = "",
		sKmatrix;
	std::vector<std::string> image_root_dirs;
	std::string recordfile;
	std::string sSfM_Data_Filename;
	std::string sMatchFilename;
	
	int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;
	std::string sIntrinsic_refinement_options = "NONE";
	bool b_use_motion_priors = false;//使用gps时记得开启
	bool b_use_pose_priors = false;//使用先验pose时开启
	EPairMode ePairmode= PAIR_FOR_MULTICAMERAS;
	//sfm_init_ImageListing	
	std::pair<bool, Vec3> gps_info;
	std::pair<bool, Pose3> pose_info;
	std::string sGPSfile;//外部GPS先验信息
	std::string sPriorWeights;	
	bool b_Group_camera_model = true;	
	int Group_camera_num = 1;
	int main_cam = 0;
	int station_num = 1;
	double focal_pixels = -1.0;	
	//Image Listing
	//Images from different lens should be stored in different folders.
	//In other words,there is a one-to-one correspondence between a len and a folder.
	int sfm_init_MCImageListing();

	//computeFeatures 
	bool bUpRight = false;
	bool bForce_f = false;
	//Describer Methods:SIFT,SIFT_ANATOMY,AKAZE_FLOAT,AKAZE_MLDB
	std::string sImage_Describer_Method = "SIFT_GPU";	
	std::string sFeaturePreset = "NORMAL";
#ifdef OPENMVG_USE_OPENMP
	int iNumThreads = 0;
#endif	 
	int computeMCFeatures();

    //computeMatches
	//Each len's K matrix is knowns
	std::string sGeometricModel = "e";
	float fDistRatio = 0.8f;
	int iMatchingVideoMode = -1;
	
	std::string sNearestMatchingMethod = "FASTCASCADEHASHINGL2";
	bool bForce_m = false;
	bool bGuided_matching = false;
	int imax_iteration = 2048;
	unsigned int ui_max_cache_size = 0;
	int computeMCMatches();

	//Pose initialization
	//globalSfM
	int iRotationAveragingMethod = int(ROTATION_AVERAGING_L2);
	int iTranslationAveragingMethod = int(TRANSLATION_AVERAGING_SOFTL1);
	vector<double> MCRSMEs;
	int globalMCSfM();

	//incrementalSfM2
	std::string sSfMInitializer_method = "STELLAR";
	int incrementalMCSfM2();

	//Eye-to-eye Calibration
	bool Evaluate_InitialPoses();
	bool Solve_AX_XB(vector<pair<int, int>>& pose_pairs, Pose3& transformation);
	bool Initial_calibration(std::vector<Pose3>& relative_poses, int main_cam, int stations);
	

	//localization
	std::string sMatchesOutDir;
	std::string sQueryDir;
	std::string sLocOutDir;
	double dMaxResidualError = 4.0;
	int Localize_camera_model = cameras::CAMERA_SPHERICAL;
	bool bUseSingleIntrinsics = true;
	bool bExportStructure = false;
	int localization();

	//Pose Refinement
	bool RBA();

	//Entrance
	int groupSfM();
};