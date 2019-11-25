
#include"fatCSFM.hpp"

bool checkGroupIntrinsicStringValidity(const std::string & Kmatrix, std::vector<double> & focal, std::vector<double> & ppx, std::vector<double> & ppy)
{
	std::vector<std::string> vec_str;
	stl::split(Kmatrix, ';', vec_str);
	if (vec_str.size() != 3*focal.size()) {
		std::cout << "\n Missing ';' character" << std::endl;
		return false;
	}
	// Check that all K matrix value are valid numbers
	for (size_t i = 0; i < vec_str.size(); ++i) {
		double readvalue = 0.0;
		std::stringstream ss;
		ss.str(vec_str[i]);
		if (!(ss >> readvalue)) {
			std::cout << "\n Used an invalid not a number character" << std::endl;
			return false;
		}
		if (i%3 == 0) focal[i/3] = readvalue;
		if (i%3 == 1) ppx[i/3] = readvalue;
		if (i%3 == 2) ppy[i/3] = readvalue;
	}
	return true;
}

double Generate_SfM_RMSE(const SfM_Data & sfm_data)
{
	IndexT residualCount = 0;
	Hash_Map<IndexT, std::vector<double>> residuals_per_view;
	std::map<IndexT, IndexT> track_length_occurences;
	for (const auto & iterTracks : sfm_data.GetLandmarks())
	{
		const Observations & obs = iterTracks.second.obs;
		track_length_occurences[obs.size()] += 1;
		for (const auto & itObs : obs)
		{
			const View * view = sfm_data.GetViews().at(itObs.first).get();
			const geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
			const cameras::IntrinsicBase * intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
			// Use absolute values
			const Vec2 residual = intrinsic->residual(pose(iterTracks.second.X), itObs.second.x).array().abs();
			residuals_per_view[itObs.first].push_back(residual(0));
			residuals_per_view[itObs.first].push_back(residual(1));
			++residualCount;
		}
	}

	// combine all residual values into one vector
	{
		IndexT residualCount = 0;
		for (const auto & residual_per_view_it : residuals_per_view)
		{
			residualCount += residual_per_view_it.second.size();
		}
		// Concat per view residual values into one vector
		std::vector<double> residuals(residualCount);
		residualCount = 0;
		for (const auto & residual_per_view_it : residuals_per_view)
		{
			std::copy(residual_per_view_it.second.begin(),
				residual_per_view_it.second.end(),
				residuals.begin() + residualCount);
			residualCount += residual_per_view_it.second.size();
		}
		if (!residuals.empty())
		{
			// RMSE computation
			const Eigen::Map<Eigen::RowVectorXd> residuals_mapping(&residuals[0], residuals.size());
			const double RMSE = std::sqrt(residuals_mapping.squaredNorm() / (double)residuals.size());
			return RMSE;
		}
	}
}

//int CompareCoords_2sfm_data()
//{
//	openMVG::sfm::SfM_Data sfm_data_gt, sfm_data_to_compare;
//	string sfmCC_file = "E:\\qj\\test\\output_outdoor\\sfm_data_psqzy.json";
//	string sfmPS_file = "E:\\qj\\test\\output_outdoor\\sfm_data_out.json";
//	if (!Load(sfm_data_gt, sfmCC_file, ESfM_Data(VIEWS | INTRINSICS | EXTRINSICS))) {
//		std::cerr << std::endl
//			<< "The input SfM_Data file \"" << sfmCC_file << "\" cannot be read." << std::endl;
//		return EXIT_FAILURE;
//	}
//	if (!Load(sfm_data_to_compare, sfmPS_file, ESfM_Data(VIEWS | INTRINSICS | EXTRINSICS))) {
//		std::cerr << std::endl
//			<< "The input SfM_Data file \"" << sfmPS_file << "\" cannot be read." << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	// Collect CC pose_ids.
//	std::set<IndexT> pose_id_gt;
//	std::transform(
//		sfm_data_gt.GetPoses().cbegin(),
//		sfm_data_gt.GetPoses().cend(),
//		std::inserter(pose_id_gt, pose_id_gt.begin()),
//		stl::RetrieveKey());
//	std::cout << "CC poses:" << pose_id_gt.size() << std::endl;
//
//	// Collect PS to compare pose_ids.
//	std::set<IndexT> pose_id_to_compare;
//	std::transform(
//		sfm_data_to_compare.GetPoses().cbegin(),
//		sfm_data_to_compare.GetPoses().cend(),
//		std::inserter(pose_id_to_compare, pose_id_to_compare.begin()),
//		stl::RetrieveKey());
//	std::cout << "PS poses:" << pose_id_to_compare.size() << std::endl;
//
//	// Check if the pose_id intersect or not
//	std::vector<IndexT> pose_id_intersection;
//	std::set_intersection(pose_id_gt.cbegin(), pose_id_gt.cend(),
//		pose_id_to_compare.cbegin(), pose_id_to_compare.cend(),
//		std::back_inserter(pose_id_intersection));
//	std::cout << "pose_id_intersection.size()=" << pose_id_intersection.size();
//
//	if (pose_id_gt.empty() || pose_id_to_compare.empty()
//		|| pose_id_intersection.size() != sfm_data_gt.GetPoses().size())
//	{
//		std::cerr << "Invalid data input. "
//			<< "The dataset does not have corresponding camera Ids." << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	std::vector<Vec3> camera_pos_gt, camera_pos_to_compare;
//	std::vector<Mat3> camera_rot_gt, camera_rot_to_compare;
//
//	for (const auto & pose_gt_it : sfm_data_gt.GetPoses())
//	{
//		const IndexT pose_id = pose_gt_it.first;
//		//sfm_data_gt.GetPoses().count(pose_id)
//		const auto & pose_to_compare_it = sfm_data_to_compare.GetPoses().at(pose_id);
//
//		camera_pos_gt.push_back(pose_gt_it.second.center());
//		camera_rot_gt.push_back(pose_gt_it.second.rotation());
//
//		camera_pos_to_compare.push_back(pose_to_compare_it.center());
//		camera_rot_to_compare.push_back(pose_to_compare_it.rotation());
//	}
//
//	openMVG::geometry::Similarity3 sim;
//	// Compute the registration:
//	if (camera_pos_gt.size() > 3)
//	{
//		const Mat X_SfM_Mat = Eigen::Map<Mat>(camera_pos_gt[0].data(), 3, camera_pos_gt.size());
//		const Mat X_GPS_Mat = Eigen::Map<Mat>(camera_pos_to_compare[0].data(), 3, camera_pos_to_compare.size());
//		std::cout << "X_SfM_Mat.size()=" << X_SfM_Mat.size() << std::endl;
//		std::cout << "X_GPS_Mat.size()=" << X_GPS_Mat.size() << std::endl;
//		/*geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
//		const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
//		if (lmeds_median != std::numeric_limits<double>::max())
//		{*/ // Compute the median residual error once the registration is applied
//		Vec3 t;
//		Mat3 R;
//		double S;
//		if (openMVG::geometry::FindRTS(X_GPS_Mat, X_SfM_Mat, &S, &t, &R))
//		{
//			openMVG::geometry::Refine_RTS(X_GPS_Mat, X_SfM_Mat, &S, &t, &R);
//			std::cout << "Found transform:\n"
//				<< " scale: " << S << "\n"
//				<< " rotation:\n" << R << "\n"
//				<< " translation: " << t.transpose() << std::endl;
//
//			const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t / S), S);
//			for (Vec3 & pos : camera_pos_to_compare) // Transform SfM poses for residual computation
//			{
//				pos = sim(pos);
//			}
//			Vec residual = (Eigen::Map<Mat3X>(camera_pos_gt[0].data(), 3, camera_pos_gt.size()) - Eigen::Map<Mat3X>(camera_pos_to_compare[0].data(), 3, camera_pos_to_compare.size())).colwise().norm();
//			std::sort(residual.data(), residual.data() + residual.size());
//			double pose_center_robust_fitting_error = residual(residual.size() / 2);
//			std::cout << "pose_center_robust_fitting_error=" << pose_center_robust_fitting_error << endl;
//		}
//	}
//
//	// -a. distance between camera center
//	std::vector<double> vec_residualErrorsX;
//	{
//		for (size_t i = 0; i < camera_pos_gt.size(); ++i) {
//			const double dResidualX = abs(camera_pos_gt[i][0] - camera_pos_to_compare[i][0]);
//			vec_residualErrorsX.push_back(dResidualX);
//		}
//	}
//	std::vector<double> vec_residualErrorsY;
//	{
//		for (size_t i = 0; i < camera_pos_gt.size(); ++i) {
//			const double dResidualY = abs(camera_pos_gt[i][1] - camera_pos_to_compare[i][1]);
//			vec_residualErrorsY.push_back(dResidualY);
//		}
//	}
//	std::vector<double> vec_residualErrorsZ;
//	{
//		for (size_t i = 0; i < camera_pos_gt.size(); ++i) {
//			const double dResidualZ = abs(camera_pos_gt[i][2] - camera_pos_to_compare[i][2]);
//			vec_residualErrorsZ.push_back(dResidualZ);
//		}
//	}
//
//	// -b. angle between rotation matrix
//	std::vector<double> vec_angularErrors;
//	{
//		std::vector<Mat3>::const_iterator iter1 = camera_rot_gt.begin();
//		for (std::vector<Mat3>::const_iterator iter2 = camera_rot_to_compare.begin();
//			iter2 != camera_rot_to_compare.end(); ++iter2, ++iter1) {
//			const Mat3 R1 = *iter1; //GT
//			const Mat3 R2T = *iter2; // Computed
//
//			const double angularErrorDegree = R2D(getRotationMagnitude(R1 * R2T.transpose()));
//			vec_angularErrors.push_back(angularErrorDegree);
//		}
//	}
//	std::cout << std::endl << "\nBaselineX error statistics : \n ";
//	minMaxMeanMedian<double>(vec_residualErrorsX.begin(), vec_residualErrorsX.end());
//
//	std::cout << std::endl << "\nBaselineY error statistics : \n ";
//	minMaxMeanMedian<double>(vec_residualErrorsY.begin(), vec_residualErrorsY.end());
//
//	std::cout << std::endl << "\nBaselineZ error statistics : \n ";
//	minMaxMeanMedian<double>(vec_residualErrorsZ.begin(), vec_residualErrorsZ.end());
//
//	std::cout << std::endl << "\nAngular error statistics : \n ";
//	minMaxMeanMedian<double>(vec_angularErrors.begin(), vec_angularErrors.end());
//
//	return getchar();
//}

features::EDESCRIBER_PRESET stringToEnum(const std::string & sPreset)
{
	features::EDESCRIBER_PRESET preset;
	if (sPreset == "NORMAL")
		preset = features::NORMAL_PRESET;
	else
		if (sPreset == "HIGH")
			preset = features::HIGH_PRESET;
		else
			if (sPreset == "ULTRA")
				preset = features::ULTRA_PRESET;
			else
				preset = features::EDESCRIBER_PRESET(-1);
	return preset;
}
//
enum EGeometricModel
{
	FUNDAMENTAL_MATRIX = 0,
	ESSENTIAL_MATRIX = 1,
	HOMOGRAPHY_MATRIX = 2,
	ESSENTIAL_MATRIX_ANGULAR = 3,
	ESSENTIAL_MATRIX_ORTHO = 4
};

// From 2 given image file-names, find the two corresponding index in the View list
bool computeIndexFromImageNames(
	const SfM_Data & sfm_data,
	const std::pair<std::string, std::string>& initialPairName,
	Pair& initialPairIndex)
{
	if (initialPairName.first == initialPairName.second)
	{
		std::cout << "\nInvalid image names. You cannot use the same image to initialize a pair." << std::endl;
		return false;
	}

	initialPairIndex = { UndefinedIndexT, UndefinedIndexT };

	/// List views filenames and find the one that correspond to the user ones:
	for (Views::const_iterator it = sfm_data.GetViews().begin();
		it != sfm_data.GetViews().end(); ++it)
	{
		const View * v = it->second.get();
		const std::string filename = stlplus::filename_part(v->s_Img_path);
		if (filename == initialPairName.first)
		{
			initialPairIndex.first = v->id_view;
		}
		else {
			if (filename == initialPairName.second)
			{
				initialPairIndex.second = v->id_view;
			}
		}
	}
	return (initialPairIndex.first != UndefinedIndexT &&
		initialPairIndex.second != UndefinedIndexT);
}

CSFM::CSFM()
{
	//std::ofstream out(recordfile);
}

CSFM::~CSFM()
{
}

enum class ESfMSceneInitializer
{
	INITIALIZE_EXISTING_POSES,
	INITIALIZE_MAX_PAIR,
	INITIALIZE_AUTO_PAIR,
	INITIALIZE_STELLAR
};

bool StringToEnum_ESfMSceneInitializer
(
	const std::string & str,
	ESfMSceneInitializer & scene_initializer
)
{
	const std::map<std::string, ESfMSceneInitializer> string_to_enum_mapping =
	{
	  {"EXISTING_POSE", ESfMSceneInitializer::INITIALIZE_EXISTING_POSES},
	  {"MAX_PAIR", ESfMSceneInitializer::INITIALIZE_MAX_PAIR},
	  {"AUTO_PAIR", ESfMSceneInitializer::INITIALIZE_AUTO_PAIR},
	  {"STELLAR", ESfMSceneInitializer::INITIALIZE_STELLAR},
	};
	auto it = string_to_enum_mapping.find(str);
	if (it == string_to_enum_mapping.end())
		return false;
	scene_initializer = it->second;
	return true;
}

int CSFM::sfm_init_MCImageListing()
{
	std::cout << " You called : " << std::endl
		<< "sfm_init_GroupImageListing " << std::endl
		<< "--imageDirectory " << sImageDir << std::endl
		<< "--outputDirectory " << sOutputDir << std::endl
		<< "--focal " << focal_pixels << std::endl
		<< "--intrinsics " << sKmatrix << std::endl
		<< "--camera_model " << i_User_camera_model << std::endl
		<< "--group_camera_num" << Group_camera_num << std::endl;

	std::pair<bool, Vec3> prior_w_info(false, Vec3(1.0, 1.0, 1.0));
	double width = -1, height = -1;
	std::vector<double> focal(Group_camera_num, -1), ppx(Group_camera_num, -1), ppy(Group_camera_num, -1);
	const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

	if (!stlplus::folder_exists(sImageDir))
	{
		std::cout << "\nThe input directory doesn't exist" << std::endl;
		return EXIT_FAILURE;
	}

	if (sOutputDir.empty())
	{
		std::cout << "\nInvalid output directory" << std::endl;
		return EXIT_FAILURE;
	}

	if (!stlplus::folder_exists(sOutputDir))
	{
		if (!stlplus::folder_create(sOutputDir))
		{
			std::cout << "\nCannot create output directory" << std::endl;
			return EXIT_FAILURE;
		}
	}

	if (sKmatrix.size() > 0 && !checkGroupIntrinsicStringValidity(sKmatrix, focal, ppx, ppy))
	{
		std::cout << "\nInvalid K matrix input" << std::endl;
		return EXIT_FAILURE;
	}
	else
		std::cout << "Kmatrix size:" << focal.size() << std::endl;

	//存储影像根目录
	sfm_data.s_root_path = sImageDir; // Setup main image root_path
	image_root_dirs = stlplus::folder_subdirectories(sImageDir);
	int i_key = 0;
	std::cout << "Total sub_directories num of sImageDir:" << image_root_dirs.size() << std::endl;
	if (image_root_dirs.size() < 1)
		return EXIT_FAILURE;

	//影像按镜头分组
	for (int subdir = 0; subdir < image_root_dirs.size(); subdir++, i_key++)
	{
		std::vector<std::string> vec_image = stlplus::folder_files(sfm_data.s_root_path + "/" + image_root_dirs[subdir]);
		std::sort(vec_image.begin(), vec_image.end());
		// Configure an empty scene with Views and their corresponding cameras

		Views & views = sfm_data.views;
		Intrinsics & intrinsics = sfm_data.intrinsics;

		std::ostringstream error_report_stream;
		for (std::vector<std::string>::const_iterator iter_image = vec_image.begin();
			iter_image != vec_image.end();
			++iter_image)
		{
			const std::string sImageFilename = stlplus::create_filespec(sImageDir + "/" + image_root_dirs[subdir], *iter_image);
			const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

			// Test if the image format is supported:
			if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown)
			{
				error_report_stream
					<< sImFilenamePart << ": Unkown image file format." << "\n";
				continue; // image cannot be opened
			}
			if (sImFilenamePart.find("mask.png") != std::string::npos
				|| sImFilenamePart.find("_mask.png") != std::string::npos)
			{
				error_report_stream
					<< sImFilenamePart << " is a mask image" << "\n";
				continue;
			}

			ImageHeader imgHeader;
			if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
				continue; // image cannot be read

			width = imgHeader.width;
			height = imgHeader.height;

			// Build intrinsic parameter related to the view
			std::shared_ptr<IntrinsicBase> intrinsic;

			if (focal[i_key] > 0 && ppx[i_key] > 0 && ppy[i_key] > 0 && width > 0 && height > 0)
			{
				// Create the desired camera type
				switch (e_User_camera_model)
				{
				case PINHOLE_CAMERA:
					intrinsic = std::make_shared<Pinhole_Intrinsic>
						(width, height, focal[i_key], ppx[i_key], ppy[i_key]);
					break;
				case PINHOLE_CAMERA_RADIAL1:
					intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
						(width, height, focal[i_key], ppx[i_key], ppy[i_key], 0.0); // setup no distortion as initial guess
					break;
				case PINHOLE_CAMERA_RADIAL3:
					intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
						(width, height, focal[i_key], ppx[i_key], ppy[i_key], 0.0, 0.0, 0.0);  // setup no distortion as initial guess
					break;
				case PINHOLE_CAMERA_BROWN:
					intrinsic = std::make_shared<Pinhole_Intrinsic_Brown_T2>
						(width, height, focal[i_key], ppx[i_key], ppy[i_key], 0.0, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
					break;
				case PINHOLE_CAMERA_FISHEYE:
					intrinsic = std::make_shared<Pinhole_Intrinsic_Fisheye>
						(width, height, focal[i_key], ppx[i_key], ppy[i_key], 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
					break;
				default:
					std::cout << "Error: unknown camera model: " << (int)e_User_camera_model << std::endl;
					return EXIT_FAILURE;
				}
			}

			View v(*iter_image, views.size(), views.size(), views.size(), width, height);

			// Add intrinsic related to the image (if any)
			if (intrinsic == nullptr)
			{
				//Since the view have invalid intrinsic data
				// (export the view, with an invalid intrinsic field value)
				v.id_intrinsic = UndefinedIndexT;
			}
			else
			{
				// Add the defined intrinsic to the sfm_container
				intrinsics[v.id_intrinsic] = intrinsic;
			}

			// Add the view to the sfm_container
			views[v.id_view] = std::make_shared<View>(v);
		}
		// Display saved warning & error messages if any.
		if (!error_report_stream.str().empty())
		{
			std::cout
				<< "\nWarning & Error messages:" << std::endl
				<< error_report_stream.str() << std::endl;
		}
	}

	if (i_key != Group_camera_num)
	{
		std::cout << "Error: num of cameras in group:" << i_key << " is not equal to Group_camera_num : " << Group_camera_num << std::endl;
		return EXIT_FAILURE;
	}
	// Group camera that share common properties if desired (leads to more faster & stable BA).
	if (b_Group_camera_model)
	{
		GroupSharedIntrinsics(sfm_data);
	}

	// Store SfM_Data views & intrinsic data
	sSfM_Data_Filename = sOutputDir + "/sfm_data.json";
	if (!Save(
		sfm_data,
		sSfM_Data_Filename,
		ESfM_Data(VIEWS | INTRINSICS)))
	{
		return EXIT_FAILURE;
	}
	station_num = sfm_data.GetViews().size() / Group_camera_num;
	std::cout << std::endl
		<< "SfMInit_ImageListing report:\n"
		<< "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
		<< "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size() << std::endl;
	return EXIT_SUCCESS;
}

int CSFM::computeMCFeatures()
{
	std::ofstream out(recordfile,ios::app);
	
	std::cout << " You called : " << std::endl
		<< "computeFeatures" << std::endl
		<< "--input_file " << sSfM_Data_Filename << std::endl
		<< "--outdir " << sOutputDir << std::endl
		<< "--describerMethod " << sImage_Describer_Method << std::endl
		<< "--upright " << bUpRight << std::endl
		<< "--describerPreset " << (sFeaturePreset.empty() ? "NORMAL" : sFeaturePreset) << std::endl
		<< "--force " << bForce_f << std::endl
		<< "--Group_camera_nums" << Group_camera_num << std::endl
#ifdef OPENMVG_USE_OPENMP
		<< "--numThreads " << iNumThreads << std::endl
#endif
		<< std::endl;
	//---------------------------------------
  // a. Load input scene
  //---------------------------------------
	//SfM_Data sfm_data;
	if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS | INTRINSICS))) {
		std::cout << std::endl
			<< "The input file \"" << sSfM_Data_Filename << "\" cannot be read" << std::endl;
		return EXIT_FAILURE;
	}

	// b. Init the image_describer
	// - retrieve the used one in case of pre-computed features
	// - else create the desired one

	using namespace openMVG::features;
	std::unique_ptr<Image_describer> image_describer;

	const std::string sImage_describer = stlplus::create_filespec(sOutputDir, "image_describer", "json");
	if (!bForce_f && stlplus::is_file(sImage_describer))
	{
		// Dynamically load the image_describer from the file (will restore old used settings)
		std::ifstream stream(sImage_describer.c_str());
		if (!stream.is_open())
			return EXIT_FAILURE;

		try
		{
			cereal::JSONInputArchive archive(stream);
			archive(cereal::make_nvp("image_describer", image_describer));
		}
		catch (const cereal::Exception & e)
		{
			std::cout << e.what() << std::endl
				<< "Cannot dynamically allocate the Image_describer interface." << std::endl;
			return EXIT_FAILURE;
		}
	}
	else
	{
		// Create the desired Image_describer method.
		// Don't use a factory, perform direct allocation
		if (sImage_Describer_Method == "SIFT")
		{
			image_describer.reset(new SIFT_Image_describer
			(SIFT_Image_describer::Params(), !bUpRight));
		}
		else if (sImage_Describer_Method == "SIFT_GPU")
			{
				image_describer.reset(new SIFT_GPU_Image_describer());
			}
		else
			if (sImage_Describer_Method == "SIFT_ANATOMY")
			{
				image_describer.reset(
					new SIFT_Anatomy_Image_describer(SIFT_Anatomy_Image_describer::Params()));
			}
			else
				if (sImage_Describer_Method == "AKAZE_FLOAT")
				{
					image_describer = AKAZE_Image_describer::create
					(AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MSURF), !bUpRight);
				}
				else
					if (sImage_Describer_Method == "AKAZE_MLDB")
					{
						image_describer = AKAZE_Image_describer::create
						(AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MLDB), !bUpRight);
					}
		if (!image_describer)
		{
			std::cout << "Cannot create the designed Image_describer:"
				<< sImage_Describer_Method << "." << std::endl;
			return EXIT_FAILURE;
		}
		else
		{
			if (!sFeaturePreset.empty())
				if (!image_describer->Set_configuration_preset(stringToEnum(sFeaturePreset)))
				{
					std::cout << "Preset configuration failed." << std::endl;
					return EXIT_FAILURE;
				}
		}

		// Export the used Image_describer and region type for:
		// - dynamic future regions computation and/or loading
		{
			std::ofstream stream(sImage_describer.c_str());
			if (!stream.is_open())
				return EXIT_FAILURE;

			cereal::JSONOutputArchive archive(stream);
			archive(cereal::make_nvp("image_describer", image_describer));
			auto regionsType = image_describer->Allocate();
			archive(cereal::make_nvp("regions_type", regionsType));
		}
	}

	// Feature extraction routines
	// For each View of the SfM_Data container:
	// - if regions file exists continue,
	// - if no file, compute features
	{
		system::Timer timer;
		Image<unsigned char> imageGray;

		C_Progress_display my_progress_bar(sfm_data.GetViews().size(),
			std::cout, "\n- EXTRACT FEATURES -\n");

		// Use a boolean to track if we must stop feature extraction
		std::atomic<bool> preemptive_exit(false);
#ifdef OPENMVG_USE_OPENMP
		const int nb_max_thread = omp_get_max_threads();

		if (iNumThreads > 0) {
			omp_set_num_threads(iNumThreads);
		}
		else {
			omp_set_num_threads(nb_max_thread);
			std::cout << "使用最大线程，当前线程数目为= " << nb_max_thread << std::endl;
		}
		
#pragma omp parallel for schedule(dynamic) if (iNumThreads > 0) private(imageGray)
#endif
		for (int i = 0; i < static_cast<int>(sfm_data.views.size()); ++i)
		{
			Views::const_iterator iterViews = sfm_data.views.begin();
			std::advance(iterViews, i);
			const View * view = iterViews->second.get();
			std::string sView_filename ;
			if (image_root_dirs.size() > 1)
				sView_filename = stlplus::create_filespec(sfm_data.s_root_path + "/" + image_root_dirs[i / station_num], view->s_Img_path);
			
			const std::string
				sFeat = stlplus::create_filespec(sOutputDir, stlplus::basename_part(sView_filename), "feat"),
				sDesc = stlplus::create_filespec(sOutputDir, stlplus::basename_part(sView_filename), "desc");

			// If features or descriptors file are missing, compute them
			if (!preemptive_exit && (bForce_f || !stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc)))
			{
				if (!ReadImage(sView_filename.c_str(), &imageGray))
					continue;

				//
				// Look if there is occlusion feature mask
				//
				Image<unsigned char> * mask = nullptr; // The mask is null by default

				const std::string
					mask_filename_local =
					stlplus::create_filespec(sfm_data.s_root_path,
						stlplus::basename_part(sView_filename) + "_mask", "png"),
					mask__filename_global =
					stlplus::create_filespec(sfm_data.s_root_path, "mask", "png");

				Image<unsigned char> imageMask;
				// Try to read the local mask
				if (stlplus::file_exists(mask_filename_local))
				{
					if (!ReadImage(mask_filename_local.c_str(), &imageMask))
					{
						std::cout << "Invalid mask: " << mask_filename_local << std::endl
							<< "Stopping feature extraction." << std::endl;
						preemptive_exit = true;
						continue;
					}
					// Use the local mask only if it fits the current image size
					if (imageMask.Width() == imageGray.Width() && imageMask.Height() == imageGray.Height())
						mask = &imageMask;
				}
				else
				{
					// Try to read the global mask
					if (stlplus::file_exists(mask__filename_global))
					{
						if (!ReadImage(mask__filename_global.c_str(), &imageMask))
						{
							std::cout << "Invalid mask: " << mask__filename_global << std::endl
								<< "Stopping feature extraction." << std::endl;
							preemptive_exit = true;
							continue;
						}
						// Use the global mask only if it fits the current image size
						if (imageMask.Width() == imageGray.Width() && imageMask.Height() == imageGray.Height())
							mask = &imageMask;
					}
				}

				// Compute features and descriptors and export them to files
				auto regions = image_describer->Describe(imageGray, mask);
				if (regions && !image_describer->Save(regions.get(), sFeat, sDesc)) {
					std::cout << "Cannot save regions for images: " << sView_filename << std::endl
						<< "Stopping feature extraction." << std::endl;
					preemptive_exit = true;
					continue;

				}
			}
			++my_progress_bar;
		}
		std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
	}
	return EXIT_SUCCESS;
}

void AdjacencyMatrixToSVG
(
	const size_t NbImages,
	const Pair_Set & corresponding_indexes,
	const std::string & sOutName
)
{
	using namespace svg;
	if (!corresponding_indexes.empty())
	{
		const float scaleFactor = 5.0f;
		svgDrawer svgStream((NbImages + 3) * 5, (NbImages + 3) * 5);
		// List possible pairs
		for (size_t I = 0; I < NbImages; ++I)
		{
			for (size_t J = 0; J < NbImages; ++J)
			{
				// If the pair have matches display a blue boxes at I,J position.
				const auto iterSearch = corresponding_indexes.find(std::make_pair(I, J));
				if (iterSearch != corresponding_indexes.end())
				{
					svgStream.drawSquare(J*scaleFactor, I*scaleFactor, scaleFactor / 2.0f,
						svgStyle().fill("blue").noStroke());
				}
			}
		}
		// Display axes with 0 -> NbImages annotation : _|
		std::ostringstream osNbImages;
		osNbImages << NbImages;
		svgStream.drawText((NbImages + 1)*scaleFactor, scaleFactor, scaleFactor, "0", "black");
		svgStream.drawText((NbImages + 1)*scaleFactor,
			(NbImages)*scaleFactor - scaleFactor, scaleFactor, osNbImages.str(), "black");
		svgStream.drawLine((NbImages + 1)*scaleFactor, 2 * scaleFactor,
			(NbImages + 1)*scaleFactor, (NbImages)*scaleFactor - 2 * scaleFactor,
			svgStyle().stroke("black", 1.0));

		svgStream.drawText(scaleFactor, (NbImages + 1)*scaleFactor, scaleFactor, "0", "black");
		svgStream.drawText((NbImages)*scaleFactor - scaleFactor,
			(NbImages + 1)*scaleFactor, scaleFactor, osNbImages.str(), "black");
		svgStream.drawLine(2 * scaleFactor, (NbImages + 1)*scaleFactor,
			(NbImages)*scaleFactor - 2 * scaleFactor, (NbImages + 1)*scaleFactor,
			svgStyle().stroke("black", 1.0));

		std::ofstream svgFileStream(sOutName.c_str());
		svgFileStream << svgStream.closeSvgFile().str();
	}
}


int CSFM::computeMCMatches()
{
	std::ofstream out(recordfile, ios::app);
	std::cout << " You called : " << "\n"
		<< "openMVG_main_ComputeMatches" << "\n"
		<< "--input_file " << sSfM_Data_Filename << "\n"
		<< "--out_dir " << sOutputDir << "\n"
		<< "Optional parameters:" << "\n"
		<< "--force " << bForce_m << "\n"
		<< "--ratio " << fDistRatio << "\n"
		<< "--geometric_model " << sGeometricModel << "\n"
		<< "--video_mode_matching " << iMatchingVideoMode << "\n"
		<< "--nearest_matching_method " << sNearestMatchingMethod << "\n"
		<< "--guided_matching " << bGuided_matching << "\n"
		<< "--cache_size " << ((ui_max_cache_size == 0) ? "unlimited" : std::to_string(ui_max_cache_size)) << std::endl;
	if (ePairmode!=PAIR_FOR_MULTICAMERAS)
		ePairmode = (iMatchingVideoMode == -1) ? PAIR_EXHAUSTIVE : PAIR_CONTIGUOUS;

	EGeometricModel eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
	std::string sGeometricMatchesFilename = "";
	switch (sGeometricModel[0])
	{
	case 'f': case 'F':
		eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
		sGeometricMatchesFilename = "matches.f.bin";
		break;
	case 'e': case 'E':
		eGeometricModelToCompute = ESSENTIAL_MATRIX;
		sGeometricMatchesFilename = "matches.e.bin";
		break;
	case 'h': case 'H':
		eGeometricModelToCompute = HOMOGRAPHY_MATRIX;
		sGeometricMatchesFilename = "matches.h.bin";
		break;
	case 'a': case 'A':
		eGeometricModelToCompute = ESSENTIAL_MATRIX_ANGULAR;
		sGeometricMatchesFilename = "matches.f.bin";
		break;
	case 'o': case 'O':
		eGeometricModelToCompute = ESSENTIAL_MATRIX_ORTHO;
		sGeometricMatchesFilename = "matches.o.bin";
		break;
	default:
		std::cout << "Unknown geometric model" << std::endl;
		return EXIT_FAILURE;
	}

	// -----------------------------
  // - Load SfM_Data Views & intrinsics data
  // a. Compute putative descriptor matches
  // b. Geometric filtering of putative matches
  // + Export some statistics
  // -----------------------------

  //---------------------------------------
  // Read SfM Scene (image view & intrinsics data)
  //---------------------------------------
	//SfM_Data sfm_data;
	if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS | INTRINSICS))) {
		std::cout << std::endl
			<< "The input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read." << std::endl;
		return EXIT_FAILURE;
	}

	//---------------------------------------
	// Load SfM Scene regions
	//---------------------------------------
	// Init the regions_type from the image describer file (used for image regions extraction)
	using namespace openMVG::features;
	const std::string sImage_describer = stlplus::create_filespec(sOutputDir, "image_describer", "json");
	std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
	if (!regions_type)
	{
		std::cout << "Invalid: "
			<< sImage_describer << " regions type file." << std::endl;
		return EXIT_FAILURE;
	}

	//---------------------------------------
	// a. Compute putative descriptor matches
	//    - Descriptor matching (according user method choice)
	//    - Keep correspondences only if NearestNeighbor ratio is ok
	//---------------------------------------

	// Load the corresponding view regions
	std::shared_ptr<Regions_Provider> regions_provider;
	if (ui_max_cache_size == 0)
	{
		// Default regions provider (load & store all regions in memory)
		regions_provider = std::make_shared<Regions_Provider>();
	}
	else
	{
		// Cached regions provider (load & store regions on demand)
		regions_provider = std::make_shared<Regions_Provider_Cache>(ui_max_cache_size);
	}

	// Show the progress on the command line:
	C_Progress_display progress;

	if (!regions_provider->load(sfm_data, sOutputDir, regions_type, &progress)) {
		std::cout << std::endl << "Invalid regions." << std::endl;
		return EXIT_FAILURE;
	}

	PairWiseMatches map_PutativesMatches;

	// Build some alias from SfM_Data Views data:
	// - List views as a vector of filenames & image sizes
	std::vector<std::string> vec_fileNames;
	std::vector<std::pair<size_t, size_t>> vec_imagesSize;
	{
		vec_fileNames.reserve(sfm_data.GetViews().size());
		vec_imagesSize.reserve(sfm_data.GetViews().size());
		for (Views::const_iterator iter = sfm_data.GetViews().begin();
			iter != sfm_data.GetViews().end();
			++iter)
		{
			const View * v = iter->second.get();
			vec_fileNames.push_back(stlplus::create_filespec(sfm_data.s_root_path + "/" + image_root_dirs[v->id_view/ station_num],
				v->s_Img_path));
			vec_imagesSize.push_back(std::make_pair(v->ui_width, v->ui_height));
		}
	}

	std::cout << std::endl << " - PUTATIVE MATCHES - " << std::endl;
	
	// If the matches already exists, reload them
	if (!bForce_m
		&& (stlplus::file_exists(sOutputDir + "/matches.putative.txt")
			|| stlplus::file_exists(sOutputDir + "/matches.putative.bin"))
		)
	{
		if (!(Load(map_PutativesMatches, sOutputDir + "/matches.putative.bin") ||
			Load(map_PutativesMatches, sOutputDir + "/matches.putative.txt")))
		{
			std::cout << "Cannot load input matches file";
			return EXIT_FAILURE;
		}
		std::cout << "\t PREVIOUS RESULTS LOADED;"
			<< " #pair: " << map_PutativesMatches.size() << std::endl;
	}
	else // Compute the putative matches
	{
		std::cout << "Use: ";
		switch (ePairmode)
		{
		case PAIR_EXHAUSTIVE: std::cout << "exhaustive pairwise matching" << std::endl; break;
		case PAIR_CONTIGUOUS: std::cout << "sequence pairwise matching" << std::endl; break;
		case PAIR_FOR_MULTICAMERAS: std::cout << "pairwise matching for multi-cameras" << std::endl; break;
		}

		// Allocate the right Matcher according the Matching requested method
		std::unique_ptr<Matcher> collectionMatcher;
		if (sNearestMatchingMethod == "AUTO")
		{
			if (regions_type->IsScalar())
			{
				std::cout << "Using FAST_CASCADE_HASHING_L2 matcher" << std::endl;
				collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
			}
			else
				if (regions_type->IsBinary())
				{
					std::cout << "Using BRUTE_FORCE_HAMMING matcher" << std::endl;
					collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_HAMMING));
				}
		}
		else
			if (sNearestMatchingMethod == "BRUTEFORCEL2")
			{
				std::cout << "Using BRUTE_FORCE_L2 matcher" << std::endl;
				collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_L2));
			}
			else
				if (sNearestMatchingMethod == "BRUTEFORCEHAMMING")
				{
					std::cout << "Using BRUTE_FORCE_HAMMING matcher" << std::endl;
					collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_HAMMING));
				}
				else
					if (sNearestMatchingMethod == "ANNL2")
					{
						std::cout << "Using ANN_L2 matcher" << std::endl;
						collectionMatcher.reset(new Matcher_Regions(fDistRatio, ANN_L2));
					}
					else
						if (sNearestMatchingMethod == "CASCADEHASHINGL2")
						{
							std::cout << "Using CASCADE_HASHING_L2 matcher" << std::endl;
							collectionMatcher.reset(new Matcher_Regions(fDistRatio, CASCADE_HASHING_L2));
						}
						else
							if (sNearestMatchingMethod == "FASTCASCADEHASHINGL2")
							{
								std::cout << "Using FAST_CASCADE_HASHING_L2 matcher" << std::endl;
								collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
							}
		if (!collectionMatcher)
		{
			std::cout << "Invalid Nearest Neighbor method: " << sNearestMatchingMethod << std::endl;
			return EXIT_FAILURE;
		}
		// Perform the matching
		system::Timer timer;
		{
			// From matching mode compute the pair list that have to be matched:
			Pair_Set pairs;
			switch (ePairmode)
			{
			case PAIR_EXHAUSTIVE: pairs = exhaustivePairs(sfm_data.GetViews().size()); break;
			case PAIR_CONTIGUOUS: pairs = contiguousWithOverlap(sfm_data.GetViews().size(), iMatchingVideoMode); break;
			case PAIR_FOR_MULTICAMERAS:pairs = MultiCamerasPairs(sfm_data.GetViews().size(), Group_camera_num); break;
			}
			// Photometric matching of putative pairs
			collectionMatcher->Match(regions_provider, pairs, map_PutativesMatches, &progress);
			//---------------------------------------
			//-- Export putative matches
			//---------------------------------------
			if (!Save(map_PutativesMatches, std::string(sOutputDir + "/matches.putative.bin")))
			{
				std::cout
					<< "Cannot save computed matches in: "
					<< std::string(sOutputDir + "/matches.putative.bin");
				return EXIT_FAILURE;
			}
		}
		std::cout << "Task (Regions Matching) done in (s): " << timer.elapsed() << std::endl;
	}
	//-- export putative matches Adjacency matrix

	PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
		map_PutativesMatches,
		stlplus::create_filespec(sOutputDir, "PutativeAdjacencyMatrix", "svg"));
	//-- export view pair graph once putative graph matches have been computed
	{
		std::set<IndexT> set_ViewIds;
		std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
			std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
		graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_PutativesMatches));
		graph::exportToGraphvizData(
			stlplus::create_filespec(sOutputDir, "putative_matches"),
			putativeGraph);
	}

	//---------------------------------------
	// b. Geometric filtering of putative matches
	//    - AContrario Estimation of the desired geometric model
	//    - Use an upper bound for the a contrario estimated threshold
	//---------------------------------------
	//如果已有geometric filter后的文件，不再重新计算
	PairWiseMatches map_GeometricMatches;

	if (!bForce_m
		&& (stlplus::file_exists(sOutputDir + "/matches.e.bin")
			|| stlplus::file_exists(sOutputDir + "/matches.f.bin"))
		)
	{
		if (!(Load(map_GeometricMatches, sOutputDir + "/matches.e.bin") ||
			Load(map_GeometricMatches, sOutputDir + "/matches.f.bin")))
		{
			std::cout << "Cannot load input matches file";
			return EXIT_FAILURE;
		}
		std::cout << "\t GeometricMatches RESULTS LOADED;"
			<< " #pair: " << map_GeometricMatches.size() << std::endl;
	}
	else 
	{
		std::unique_ptr<ImageCollectionGeometricFilter> filter_ptr(
			new ImageCollectionGeometricFilter(&sfm_data, regions_provider));

		if (filter_ptr)
		{
			system::Timer timer;
			const double d_distance_ratio = 0.6;


			switch (eGeometricModelToCompute)
			{
			case HOMOGRAPHY_MATRIX:
			{
				const bool bGeometric_only_guided_matching = true;
				filter_ptr->Robust_model_estimation(
					GeometricFilter_HMatrix_AC(4.0, imax_iteration),
					map_PutativesMatches, bGuided_matching,
					bGeometric_only_guided_matching ? -1.0 : d_distance_ratio, &progress);
				map_GeometricMatches = filter_ptr->Get_geometric_matches();
			}
			break;
			case FUNDAMENTAL_MATRIX:
			{
				filter_ptr->Robust_model_estimation(
					GeometricFilter_FMatrix_AC(4.0, imax_iteration),
					map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
				map_GeometricMatches = filter_ptr->Get_geometric_matches();
			}
			break;
			case ESSENTIAL_MATRIX:
			{
				filter_ptr->Robust_model_estimation(
					GeometricFilter_EMatrix_AC(4.0, imax_iteration),
					map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
				map_GeometricMatches = filter_ptr->Get_geometric_matches();

				//-- Perform an additional check to remove pairs with poor overlap
				std::vector<PairWiseMatches::key_type> vec_toRemove;
				for (const auto & pairwisematches_it : map_GeometricMatches)
				{
					const size_t putativePhotometricCount = map_PutativesMatches.find(pairwisematches_it.first)->second.size();
					const size_t putativeGeometricCount = pairwisematches_it.second.size();
					const float ratio = putativeGeometricCount / static_cast<float>(putativePhotometricCount);
					if (putativeGeometricCount < 50 || ratio < .3f) {
						// the pair will be removed
						vec_toRemove.push_back(pairwisematches_it.first);
					}
				}
				//-- remove discarded pairs
				for (const auto & pair_to_remove_it : vec_toRemove)
				{
					map_GeometricMatches.erase(pair_to_remove_it);
				}
			}
			break;
			case ESSENTIAL_MATRIX_ANGULAR:
			{
				filter_ptr->Robust_model_estimation(
					GeometricFilter_ESphericalMatrix_AC_Angular(4.0, imax_iteration),
					map_PutativesMatches, bGuided_matching);
				map_GeometricMatches = filter_ptr->Get_geometric_matches();
			}
			break;
			case ESSENTIAL_MATRIX_ORTHO:
			{
				filter_ptr->Robust_model_estimation(
					GeometricFilter_EOMatrix_RA(2.0, imax_iteration),
					map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
				map_GeometricMatches = filter_ptr->Get_geometric_matches();
			}
			break;
			}

			//---------------------------------------
			//-- Export geometric filtered matches
			//---------------------------------------
			if (!Save(map_GeometricMatches,
				std::string(sOutputDir + "/" + sGeometricMatchesFilename)))
			{
				std::cout
					<< "Cannot save computed matches in: "
					<< std::string(sOutputDir + "/" + sGeometricMatchesFilename);
				return EXIT_FAILURE;
			}
			std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
			
	    }
		//-- export Adjacency matrix
		std::cout << "\n Export Adjacency Matrix of the pairwise's geometric matches"
			<< std::endl;
		PairWiseMatchingToAdjacencyMatrixSVG_MC(vec_fileNames.size(),Group_camera_num,
			map_GeometricMatches,
			stlplus::create_filespec(sOutputDir, "GeometricAdjacencyMatrix", "svg"));
		//-- export view pair graph once geometric filter have been done
		{
			std::set<IndexT> set_ViewIds;
			std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
				std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
			graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_GeometricMatches));
			graph::exportToGraphvizData(
				stlplus::create_filespec(sOutputDir, "geometric_matches"),
				putativeGraph);
		}
	}

	return EXIT_SUCCESS;
}

int CSFM::globalMCSfM()
{
	std::ofstream record_out(recordfile, ios::app);
	if (iRotationAveragingMethod < ROTATION_AVERAGING_L1 ||
		iRotationAveragingMethod > ROTATION_AVERAGING_L2) {
		std::cout << "\n Rotation averaging method is invalid" << std::endl;
		return EXIT_FAILURE;
	}

	// Load input SfM_Data scene
	SfM_Data sfm_data;
	if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS | INTRINSICS | EXTRINSICS))) {
		std::cout << std::endl
			<< "The input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read." << std::endl;
		return EXIT_FAILURE;
	}

	const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
		cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);
	if (intrinsic_refinement_options == static_cast<cameras::Intrinsic_Parameter_Type>(0))
	{
		std::cout << "Invalid input for Bundle Adjusment Intrinsic parameter refinement option" << std::endl;
		return EXIT_FAILURE;
	}

	if (iTranslationAveragingMethod < TRANSLATION_AVERAGING_L1 ||
		iTranslationAveragingMethod > TRANSLATION_AVERAGING_SOFTL1) {
		std::cout << "\n Translation averaging method is invalid" << std::endl;
		return EXIT_FAILURE;
	}
	// Init the regions_type from the image describer file (used for image regions extraction)
	using namespace openMVG::features;
	const std::string sImage_describer = stlplus::create_filespec(sOutputDir, "image_describer", "json");
	std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
	if (!regions_type)
	{
		std::cout << "Invalid: "
			<< sImage_describer << " regions type file." << std::endl;
		return EXIT_FAILURE;
	}

	//---------------------------------------
  // Global SfM reconstruction process for each camera
  //---------------------------------------
	station_num = sfm_data.views.size() / Group_camera_num;
	for (int cam = 0; cam < Group_camera_num; cam++)
	{
		SfM_Data sfm_data_single;
		// 1. Fill the  single sfm_data scene with view, intrinsic
		for (int i=0;i<station_num;i++)
		{
			int view_id = cam*station_num + i;
			int intrinsic_id = sfm_data.GetViews().at(view_id)->id_intrinsic;
			int pose_id= sfm_data.GetViews().at(view_id)->id_pose;
			sfm_data_single.views.insert(*(sfm_data.GetViews().find(view_id)));
			sfm_data_single.intrinsics.insert(*(sfm_data.GetIntrinsics().find(intrinsic_id)));
		}
		
		if (sfm_data_single.views.size() != station_num)
		{
			std::cout << "Not valid sfm_data part!!!" << std::endl;
			return EXIT_FAILURE;
		}

		// 2. Begin the single globalsfm
		openMVG::system::Timer timer;
		GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
			sfm_data_single,
			sOutputDir,
			stlplus::create_filespec(sOutputDir, "Reconstruction_Report.html"));

		// Features reading
		std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
		if (!feats_provider->load(sfm_data_single, sOutputDir, regions_type)) {
			std::cout << std::endl
				<< "Invalid features." << std::endl;
			return EXIT_FAILURE;
		}
		// Matches reading
		std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
		if // Try to read the provided match filename or the default one (matches.e.txt/bin)
			(
				!(matches_provider->load(sfm_data_single, sMatchFilename) ||
					matches_provider->load(sfm_data_single, stlplus::create_filespec(sOutputDir, "matches.e.bin")) ||
					matches_provider->load(sfm_data_single, stlplus::create_filespec(sOutputDir, "matches.f.bin")))
				)
		{
			std::cout << std::endl
				<< "Invalid matches file." << std::endl;
			return EXIT_FAILURE;
		}

		// Configure the features_provider & the matches_provider
		sfmEngine.SetFeaturesProvider(feats_provider.get());
		sfmEngine.SetMatchesProvider(matches_provider.get());

		// Configure reconstruction parameters
		sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);

		sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);

		// Configure motion averaging method
		sfmEngine.SetRotationAveragingMethod(
			ERotationAveragingMethod(iRotationAveragingMethod));
		sfmEngine.SetTranslationAveragingMethod(
			ETranslationAveragingMethod(iTranslationAveragingMethod));

		if (sfmEngine.Process())
		{
			MCRMSEs[cam] = Generate_SfM_RMSE(sfm_data_single);
			//-- Export to disk computed scene (data & visualizable results)
			std::cout << "...Export SfM_Data to disk." << std::endl;
			Save(sfmEngine.Get_SfM_Data(),
				stlplus::create_filespec(sOutputDir, "sfm_data_for_cam_"+to_string(cam), ".bin"),
				ESfM_Data(ALL));

			Save(sfmEngine.Get_SfM_Data(),
				stlplus::create_filespec(sOutputDir, "cloud_and_poses_for_cam_" + to_string(cam), ".ply"),
				ESfM_Data(ALL));
		}
		else
		{
			std::cout << "Failed during sfmEngine for cam"<<to_string(cam)<<" processing!!!" << std::endl;
			continue;
		}
		
		//3.Integrate the single sfm_data scene poses and structure into the global sfm_data scene poses
		for (auto & pose_it : sfmEngine.Get_SfM_Data().GetPoses())
			sfm_data.poses[pose_it.first] = pose_it.second;
		record_out << "After Cam_" << to_string(cam) << " incremental_v2 sfm, poses size:" << sfmEngine.Get_SfM_Data().poses.size() << "/" << sfm_data.poses.size() << std::endl;

		for (auto & landmark_it : sfmEngine.Get_SfM_Data().GetLandmarks())
		{
			int landmark_id = sfm_data.structure.size();
			Landmark landmark = landmark_it.second;
			sfm_data.structure[landmark_id] = landmark;
		}
		record_out << "After Cam_" << to_string(cam) << " incremental_v2 sfm, structure size:" << sfmEngine.Get_SfM_Data().structure.size() << "/" << sfm_data.structure.size() << std::endl;

	}

	Save(sfm_data,
		stlplus::create_filespec(sOutputDir, "sfm_data_intial", ".bin"),
		ESfM_Data(ALL));

	return EXIT_SUCCESS;
}

int CSFM::incrementalMCSfM2()
{
	std::ofstream record_out(recordfile, ios::app);
	std::cout << "Sequential/Incremental reconstruction (Engine 2)" << std::endl
		<< std::endl;

	const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
		cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);
	if (intrinsic_refinement_options == static_cast<cameras::Intrinsic_Parameter_Type>(0))
	{
		std::cout << "Invalid input for the Bundle Adjusment Intrinsic parameter refinement option" << std::endl;
		return EXIT_FAILURE;
	}

	ESfMSceneInitializer scene_initializer_enum;
	if (!StringToEnum_ESfMSceneInitializer(sSfMInitializer_method, scene_initializer_enum))
	{
		std::cout << "Invalid input for the SfM initializer option" << std::endl;
		return EXIT_FAILURE;
	}

	// Load input SfM_Data scene
	SfM_Data sfm_data;
	if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS | INTRINSICS | EXTRINSICS))) {
		std::cout << std::endl
			<< "The input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read." << std::endl;
		return EXIT_FAILURE;
	}

	// Init the regions_type from the image describer file (used for image regions extraction)
	using namespace openMVG::features;
	const std::string sImage_describer = stlplus::create_filespec(sOutputDir, "image_describer", "json");
	std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
	if (!regions_type)
	{
		std::cout << "Invalid: "
			<< sImage_describer << " regions type file." << std::endl;
		return EXIT_FAILURE;
	}

	
	// Sequential reconstruction process
  //---------------------------------------
	station_num = sfm_data.views.size() / Group_camera_num;
	for (int cam = 0; cam < Group_camera_num; cam++)
	{
		SfM_Data sfm_data_single;
		// 1. Fill the  single sfm_data scene with view, intrinsic
		for (int i = 0; i < station_num; i++)
		{
			int view_id = cam*station_num + i;
			int intrinsic_id = sfm_data.GetViews().at(view_id)->id_intrinsic;
			int pose_id = sfm_data.GetViews().at(view_id)->id_pose;
			sfm_data_single.views.insert(*(sfm_data.GetViews().find(view_id)));
			sfm_data_single.intrinsics.insert(*(sfm_data.GetIntrinsics().find(intrinsic_id)));
		}

		if (sfm_data_single.views.size() != station_num)
		{
			std::cout << "Not valid sfm_data part!!!" << std::endl;
			return EXIT_FAILURE;
		}
		// Features reading
		std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
		if (!feats_provider->load(sfm_data_single, sOutputDir, regions_type)) {
			std::cout << std::endl
				<< "Invalid features." << std::endl;
			return EXIT_FAILURE;
		}
		// Matches reading
		std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
		if // Try to read the provided match filename or the default one (matches.e.txt/bin)
			(
				!(matches_provider->load(sfm_data_single, sMatchFilename) ||
					matches_provider->load(sfm_data_single, stlplus::create_filespec(sOutputDir, "matches.e.bin")) ||
					matches_provider->load(sfm_data_single, stlplus::create_filespec(sOutputDir, "matches.f.bin")))
				)
		{
			std::cout << std::endl
				<< "Invalid matches file." << std::endl;
			return EXIT_FAILURE;
		}
		openMVG::system::Timer timer;

		std::unique_ptr<SfMSceneInitializer> scene_initializer;
		switch (scene_initializer_enum)
		{
		case ESfMSceneInitializer::INITIALIZE_AUTO_PAIR:
			std::cout << "Not yet implemented." << std::endl;
			return EXIT_FAILURE;
			break;
		case ESfMSceneInitializer::INITIALIZE_MAX_PAIR:
			scene_initializer.reset(new SfMSceneInitializerMaxPair(sfm_data_single,
				feats_provider.get(),
				matches_provider.get()));
			break;
		case ESfMSceneInitializer::INITIALIZE_EXISTING_POSES:
			scene_initializer.reset(new SfMSceneInitializer(sfm_data_single,
				feats_provider.get(),
				matches_provider.get()));
			break;
		case ESfMSceneInitializer::INITIALIZE_STELLAR:
			scene_initializer.reset(new SfMSceneInitializerStellar(sfm_data_single,
				feats_provider.get(),
				matches_provider.get()));
			break;
		default:
			return EXIT_FAILURE;
		}
		if (!scene_initializer)
		{
			std::cout << "Invalid scene initializer." << std::endl;
			return EXIT_FAILURE;
		}

		SequentialSfMReconstructionEngine2 sfmEngine(
			scene_initializer.get(),
			sfm_data_single,
			sOutputDir,
			stlplus::create_filespec(sOutputDir, "Reconstruction_Report_cam" + to_string(cam)+".html"));

		// Configure the features_provider & the matches_provider
		sfmEngine.SetFeaturesProvider(feats_provider.get());
		sfmEngine.SetMatchesProvider(matches_provider.get());

		// Configure reconstruction parameters
		sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
		sfmEngine.SetUnknownCameraType(EINTRINSIC(i_User_camera_model));
		sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);

		if (sfmEngine.Process())
		{
			std::cout << std::endl << " Total Ac-Sfm took (s): " << timer.elapsed() << std::endl;
			//-- Export to disk computed scene (data & visualizable results)
			std::cout << "...Export SfM_Data to disk." << std::endl;
			Save(sfmEngine.Get_SfM_Data(),
				stlplus::create_filespec(sOutputDir, "sfm_data_for_cam_" + to_string(cam), ".json"),
				ESfM_Data(EXTRINSICS|STRUCTURE));

			Save(sfmEngine.Get_SfM_Data(),
				stlplus::create_filespec(sOutputDir, "cloud_and_poses_for_cam_" + to_string(cam), ".ply"),
				ESfM_Data(ALL));
		}
		else
		{
			std::cout << "Failed during sfmEngine for cam" << to_string(cam) << " processing!!!" << std::endl;
			continue;
		}

		//3.Integrate the single sfm_data scene poses and structure into the global sfm_data scene poses
		for (auto & pose_it : sfmEngine.Get_SfM_Data().GetPoses())
			sfm_data.poses[pose_it.first] = pose_it.second;
		record_out<< "After Cam_" << to_string(cam) << " incremental_v2 sfm, poses size:" << sfmEngine.Get_SfM_Data().poses.size()<<"/"<<sfm_data.poses.size() << std::endl;

		for (auto & landmark_it : sfmEngine.Get_SfM_Data().GetLandmarks())
		{
			int landmark_id = sfm_data.structure.size();
			Landmark landmark = landmark_it.second;
			sfm_data.structure[landmark_id] = landmark;
		}
		record_out << "After Cam_" << to_string(cam) << " incremental_v2 sfm, structure size:" << sfmEngine.Get_SfM_Data().structure.size() << "/"<<sfm_data.structure.size() << std::endl;
	}

	Save(sfm_data,
		stlplus::create_filespec(sOutputDir, "sfm_data_intial", ".bin"),
		ESfM_Data(ALL));
	return EXIT_SUCCESS;
}

//many pairs of observations {(A1,B1),(A2,B2),⋯,(Ak,Bk)}，finding the soluion of the eqution AX=XB can be turned into minimization problem
bool CSFM::Solve_AX_XB(vector<pair<int, int>>& pose_pairs, Pose3& transformation)
{
	if (pose_pairs.size() < 2)
	{
		std::cout << "Not enough pose pairs!!!!!" << std::endl;
		return false;
	}
	return true;
}

bool CSFM::Evaluate_InitialPoses()
{
	string sfm_Data_Filename = stlplus::create_filespec(sOutputDir, "sfm_data_intial", ".bin");
	if (!Load(sfm_data, sfm_Data_Filename, ESfM_Data(ALL))) {
		std::cout << std::endl
			<< "The input SfM_Data file \"" << sfm_Data_Filename << "\" cannot be read." << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << std::endl
		<< "The input SfM_Data file \n" 
		<< " #views: " << sfm_data.views.size() << "\n"
		<< " #poses: " << sfm_data.poses.size() << "\n"
		<< " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
		<< " #tracks: " << sfm_data.structure.size() << "\n"
		<< std::endl;

	//Count number of poses for each camera
	vector<int> posesnum_each_camera;//Poses' num of each camera's sfm_data
	posesnum_each_camera.resize(Group_camera_num, 0);

	//Build the station information matrix which can be used to judge if the poses for each camera are valid in the same station.
	std::vector<std::vector<bool>> station_info ;//Station matrix information
	for (int i = 0; i < station_num; i++)
	{
		std::vector<bool> valid_cams;
		valid_cams.resize(Group_camera_num, 0);
		for (int cam = 0; cam < Group_camera_num; cam++)
		{
			int view_id = cam*station_num + i;
			const View* view = sfm_data.GetViews().at(view_id).get();
			if (sfm_data.IsPoseAndIntrinsicDefined(view))
			{
				valid_cams[cam] = 1;
				posesnum_each_camera[cam] += 1;
			}
		}
		station_info.push_back(valid_cams);
	}

	//Build pairs for solve  AX=BX problems.
	map<pair<int,int>,vector<pair<int,int>>>main_and_slaves;//(main_cam_id,slave_cam_id)-->[main_cam_poseid,slave_cam_poseid]
	map<pair<int, int>, Pose3> transforms;
	for (int main_it = 0; main_it < Group_camera_num; main_it++)
	{
		for (int slave_it = main_it + 1; slave_it < Group_camera_num; slave_it++)
		{
			vector<pair<int, int>> pose_pairs;
			for (int i = 0; i < station_num; i++)
			{
				if (station_info[i][main_it] && station_info[i][slave_it])
					pose_pairs.push_back(make_pair(main_it*station_num + i, slave_it*station_num + i));
			}
			main_and_slaves[make_pair(main_it, slave_it)] = pose_pairs;
			pose_pairs.clear();
		}
	}

	//solve
	for(auto & transform_it:main_and_slaves)
	{
		//At least 2 pairs
		Pose3 transformation = Pose3();
		if (!Solve_AX_XB(transform_it.second, transformation))
		{
			std::cout << "Unable to solve the transformation between camera_" << transform_it.first.first << "and camera_" << transform_it.first.second << std::endl;
			continue;
		}
		transforms[transform_it.first] = transformation;
		std::cout << "Able to solve the transformation between camera_" << transform_it.first.first << "and camera_" << transform_it.first.second << std::endl;
	}
	
	//Assert the main cam with minimal sfm_data RMSE
	for (int cam = 0; cam < Group_camera_num; cam++)
		if (MCRMSEs[cam] < MCRMSEs[main_cam])
			main_cam = cam;
	std::cout << "We choose camera_" << main_cam << " as our main camera with minimal SfM RMSE:" << MCRMSEs[main_cam] << std::endl;
	return true;
}
//bool CSFM::Initial_calibration(std::vector<Pose3>& relative_poses,int main_cam,int stations)
//{
//	//用平均偏移量作为标定值
//			std::vector<Vec3> cen_sum;
//			cen_sum.resize(Group_camera_num, Vec3(0.0, 0.0, 0.0));
//			std::vector<Vec4> quat_sum;
//			quat_sum.resize(Group_camera_num, Vec4(0.0, 0.0, 0.0, 0.0));
//			std::map<openMVG::IndexT, std::vector<Pose3>> relative_poses_groups;
//			for (int j= 0; j < stations; j++)
//				for (int i = 0; i < Group_camera_num; i++)
//					relative_poses_groups[j].push_back(Pose3());
//
//			for (const auto & view_it:sfm_data.views)
//			{
//				int cam_id = 0 ,station = 0, main_pose_id = 0;
//				if (b_views_grouped)
//				{
//					cam_id = (int)view_it.first%Group_camera_num;
//					station = (int)view_it.first / Group_camera_num;
//					main_pose_id = station*Group_camera_num + main_cam;
//				}
//				else
//				{
//					cam_id = (int)view_it.first / stations;
//					station = (int)view_it.first%stations;
//					main_pose_id = main_cam*stations +station;
//				}
//				Pose3 main_cam_pose = sfm_data.GetPoses().at(main_pose_id);
//				if (cam_id==main_cam)
//					continue;
//				else
//				{
//					Pose3 pose = sfm_data.GetPoseOrDie(view_it.second.get());
//					Vec3 trans = pose.center() - main_cam_pose.center();
//					Mat3 rot = pose.rotation()*(main_cam_pose.rotation().transpose());
//					relative_poses_groups[station][cam_id]=Pose3(rot, trans);
//				}
//			}
//
//			for (int i=0;i<Group_camera_num;i++)
//			{
//				for (int j = 0; j <stations ; j++)
//				{
//					cen_sum[i] = cen_sum[i] + relative_poses_groups[j][i].center();
//					Eigen::Matrix3d rot_i = relative_poses_groups[j][i].rotation();
//					Eigen::Quaterniond quat(rot_i);
//					quat_sum[i] = quat_sum[i] + Vec4(quat.x(), quat.y(), quat.z(), quat.w());
//				}
//			}
//			std::cout << "cen_sum" <<cen_sum[main_cam] << std::endl;
//
//			relative_poses.resize(Group_camera_num, Pose3());
//			for (int i = 0; i < Group_camera_num; i++)
//			{
//				if (i==main_cam)
//					continue;
//				Vec3 relative_t_average = cen_sum[i] / relative_poses_groups.size();
//				Vec4 relative_q_average = quat_sum[i] / relative_poses_groups.size();
//				Eigen::Quaterniond quat_average(relative_q_average);
//				Mat3 relative_r_average = quat_average.toRotationMatrix();
//				relative_poses[i] = Pose3(relative_r_average, relative_t_average);
//				const double angularErrorDegree = R2D(getRotationMagnitude(relative_r_average));
//				std::cout << "Camera" << i << ": " << std::endl;
//				std::cout << "Translation:" << relative_t_average[0] << "," << relative_t_average[1] << "," << relative_t_average[2] << std::endl;
//				std::cout << "Rotation:" << relative_r_average(0, 0) << "," << relative_r_average(0, 1) << "," << relative_r_average(0, 2) << ","
//					<< relative_r_average(1, 0) << "," << relative_r_average(1, 1) << "," << relative_r_average(1, 2) << ","
//					<< relative_r_average(2, 0) << "," << relative_r_average(2, 1) << "," << relative_r_average(2, 2) << std::endl;
//				std::cout << "AngularRelativeDegree:" << angularErrorDegree << std::endl;
//				std::cout << std::endl;
//			}
//
//			return relative_poses.size() == Group_camera_num;
//}

//bool CSFM::Evaluate_InitialPoses(int stations)
//{
//	//七目镜头有些影像必然会丢矢
//	//用其他全定向成功的站来推测丢失影像的位姿
//
//	//相对关系矩阵
//	//按照（i,j),i<j构建上三角矩阵存储平均值
//	std::cout << ".........................Evaluating poses from known poses..................." << std::endl;
//	std::map<std::pair<int, int>, std::vector<geometry::Pose3>> cam_cali_sum;
//	std::map<std::pair<int, int>, geometry::Pose3> cam_cali;
//
//	//遍历每一站
//	for (int s = 0; s < stations; s++)
//	{
//		for (int cam_i = 0; cam_i < Group_camera_num; cam_i++)
//			for (int cam_j = cam_i; cam_j < Group_camera_num; cam_j++)
//			{
//				if (cam_i == cam_j)
//					cam_cali_sum[std::make_pair(cam_i, cam_j)].push_back(geometry::Pose3());
//				else
//				{
//					View* view_I;
//					View* view_J;
//					if (b_views_grouped)
//					{
//						view_I = sfm_data.GetViews().at(cam_i+ s*Group_camera_num).get();
//						view_J = sfm_data.GetViews().at(cam_j +s*Group_camera_num).get();
//					}
//					else 
//					{
//						view_I = sfm_data.GetViews().at(cam_i*stations + s).get();
//						view_J = sfm_data.GetViews().at(cam_j*stations + s).get();
//					}
//
//					if (sfm_data.IsPoseAndIntrinsicDefined(view_I) && sfm_data.IsPoseAndIntrinsicDefined(view_J))
//					{
//						Pose3 pose_I = sfm_data.GetPoseOrDie(view_I);
//						Pose3 pose_J = sfm_data.GetPoseOrDie(view_J);
//						Vec3 trans = pose_J.center() - pose_I.center();
//						Mat3 rot = pose_J.rotation()*(pose_I.rotation().transpose());
//						cam_cali_sum[std::make_pair(cam_i, cam_j)].push_back(geometry::Pose3(rot, trans));
//					}
//					else
//						continue;
//				}
//			}
//	}
//	
//	//再针对每一对关系求平均
//	for (int cam_i = 0; cam_i < Group_camera_num; cam_i++)
//		for (int cam_j = cam_i; cam_j < Group_camera_num; cam_j++)
//		{
//			if (cam_i == cam_j)
//				cam_cali[std::make_pair(cam_i, cam_j)] = Pose3();
//			else
//			{
//				//如果任何一个相对关系为空则返回false
//				if (cam_cali_sum[std::make_pair(cam_i, cam_j)].empty())
//					return false;
//
//				//求poses的平均值
//				Vec3 cen_sum=Vec3(0,0,0);
//				Vec4 quat_sum= Vec4(0.0, 0.0, 0.0, 0.0);
//				std::cout << "cam_cali_sum[(" << cam_i << "," << cam_j << ")].size()=" << cam_cali_sum[std::make_pair(cam_i, cam_j)].size() << std::endl;
//				for (int count = 0; count < cam_cali_sum[std::make_pair(cam_i, cam_j)].size(); count++)
//				{
//					cen_sum+= cam_cali_sum[std::make_pair(cam_i, cam_j)][count].center();
//					Eigen::Matrix3d rot_i = cam_cali_sum[std::make_pair(cam_i, cam_j)][count].rotation();
//					Eigen::Quaterniond quat(rot_i);
//					quat_sum+= Vec4(quat.x(), quat.y(), quat.z(), quat.w());
//				}
//				Vec3 relative_t_average = cen_sum / cam_cali_sum[std::make_pair(cam_i, cam_j)].size();
//				Vec4 relative_q_average = quat_sum / cam_cali_sum[std::make_pair(cam_i, cam_j)].size();
//				Eigen::Quaterniond quat_average(relative_q_average);
//				Mat3 relative_r_average = quat_average.toRotationMatrix();
//				cam_cali[std::make_pair(cam_i, cam_j)] = Pose3(relative_r_average, relative_t_average);
//			}
//		}
//	
//	//然后遍历没有pose的view，估计它们的位姿
//	std::cout << "Get views without poses and give them evaluated average values............" <<std::endl;
//	for (const auto & view_it : sfm_data.views)
//	{
//		if (sfm_data.IsPoseAndIntrinsicDefined(view_it.second.get()))
//			continue;
//		int cam_id = 0;
//		int cam_station = 0;
//		std::vector<std::pair<int,int>> near_camid_views;
//		if (b_views_grouped)
//		{
//			cam_id = (int)view_it.first%Group_camera_num;
//			cam_station = (int)view_it.first / Group_camera_num;
//			for (int cam_i = 0; cam_i < Group_camera_num; cam_i++)
//				near_camid_views.push_back(std::make_pair(cam_i, cam_station*Group_camera_num + cam_i));
//		}
//		else
//		{
//			cam_id = view_it.first / stations;
//			cam_station = view_it.first%stations; 
//			for (int cam_i = 0; cam_i < Group_camera_num; cam_i++)
//				near_camid_views.push_back(std::make_pair(cam_i, cam_i*stations + cam_station));
//		}
//		std::cout << "Cam_id to evaluating:" << cam_id << " and near_camid_views.size()=" << near_camid_views.size() << std::endl;
//		std::vector<Pose3> evaluated_poses;
//		
//		for (int i=0;i<near_camid_views.size();i++)
//		{
//			View* v = sfm_data.views.at(near_camid_views[i].second).get();
//			if (!sfm_data.IsPoseAndIntrinsicDefined(v))
//			{
//				std::cout << "near views :" << v->id_view << " also doesn't have a difined pose!!" << std::endl;
//				continue;
//			}
//			
//			if (near_camid_views[i].first <= cam_id)
//			{
//				Pose3 cam_i_pose = sfm_data.GetPoseOrDie(v);
//				std::cout << "pose_id:" << v->id_pose << ",and pose.center is" << cam_i_pose.center()[0] << " " << cam_i_pose.center()[1] << " " << cam_i_pose.center()[2] << endl;
//				Pose3 relative_pose = cam_cali[std::make_pair(near_camid_views[i].first, cam_id)];
//				Pose3 cam_id_pose = Pose3(relative_pose.rotation()*cam_i_pose.rotation(), cam_i_pose.center() + relative_pose.center());
//				evaluated_poses.push_back(cam_id_pose);
//			}
//			else
//			{
//				Pose3 cam_i_pose = sfm_data.GetPoseOrDie(v);
//				std::cout << "pose_id:" << v->id_pose << ",and pose.center is" << cam_i_pose.center()[0] << " " << cam_i_pose.center()[1] << " " << cam_i_pose.center()[2] << endl;
//				Pose3 relative_pose = cam_cali[std::make_pair(cam_id, near_camid_views[i].first)];
//				Pose3 cam_id_pose = Pose3((relative_pose.rotation().transpose())*cam_i_pose.rotation(), cam_i_pose.center() - relative_pose.center());
//				evaluated_poses.push_back(cam_id_pose);
//			}
//		}
//		if (evaluated_poses.size() < 1)
//			return false;
//		//取平均值
//		//求poses的平均值
//		Vec3 pcen_sum = Vec3(0, 0, 0);
//		Vec4 pquat_sum = Vec4(0.0, 0.0, 0.0, 0.0);
//		for (int i = 0; i <evaluated_poses.size() ; i++)
//		{
//			pcen_sum += evaluated_poses[i].center();
//			Eigen::Matrix3d prot_i = evaluated_poses[i].rotation();
//			Eigen::Quaterniond pquat(prot_i);
//			pquat_sum += Vec4(pquat.x(), pquat.y(), pquat.z(), pquat.w());
//		}
//		Vec3 prelative_t_average = pcen_sum / evaluated_poses.size();
//		Vec4 prelative_q_average = pquat_sum / evaluated_poses.size();
//		Eigen::Quaterniond pquat_average(prelative_q_average);
//		Mat3 prelative_r_average = pquat_average.toRotationMatrix();
//		sfm_data.poses[view_it.second->id_pose]= Pose3(prelative_r_average, prelative_t_average);
//		std::cout << "Evaluated pose_id:" << view_it.second->id_pose << ",and the pose.center() is" << prelative_t_average << endl;
//	}
//	return sfm_data.views.size()==sfm_data.poses.size();
//}

//int CSFM::localization()
//{
//	std::ofstream out(recordfile, ios::ate);
//	// Load input SfM_Data scene
//	SfM_Data sfm_data;
//	if (!Load(sfm_data, sOutputDir+"/sfm_data.bin", ESfM_Data(ALL))) {
//		std::cerr << std::endl
//			<< "The input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read." << std::endl;
//		return EXIT_FAILURE;
//	}
//	out << ";Load Sfm_data #views:" << sfm_data.views.size() << ",#poses:" << sfm_data.poses.size() << std::endl;
//	out << ";Result for " << sQueryDir << std::endl;
//	out << ";folder,x_coordinates,y_coordinates,z_coordinates" << std::endl;
//
//	if (sMatchesOutDir.empty())
//	{
//		sMatchesOutDir = sLocOutDir;
//	}
//
//	if (sfm_data.GetPoses().empty() || sfm_data.GetLandmarks().empty())
//	{
//		std::cerr << std::endl
//			<< "The input SfM_Data file have not 3D content to match with." << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	// ---------------
//	// Initialization
//	// ---------------
//
//	// Init the regions_type from the image describer file (used for image regions extraction)
//	using namespace openMVG::features;
//	const std::string sImage_describer = stlplus::create_filespec(sOutputDir, "image_describer", "json");
//	std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
//	if (!regions_type)
//	{
//		std::cerr << "Invalid: "
//			<< sImage_describer << " regions type file." << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	// Init the feature extractor that have been used for the reconstruction
//	std::unique_ptr<Image_describer> image_describer;
//	if (stlplus::is_file(sImage_describer))
//	{
//		// Dynamically load the image_describer from the file (will restore old used settings)
//		std::ifstream stream(sImage_describer.c_str());
//		if (!stream.is_open())
//			return EXIT_FAILURE;
//
//		try
//		{
//			cereal::JSONInputArchive archive(stream);
//			archive(cereal::make_nvp("image_describer", image_describer));
//		}
//		catch (const cereal::Exception & e)
//		{
//			std::cerr << e.what() << std::endl
//				<< "Cannot dynamically allocate the Image_describer interface." << std::endl;
//			return EXIT_FAILURE;
//		}
//	}
//	else
//	{
//		std::cerr << "Expected file image_describer.json cannot be opened." << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	// Show the progress on the command line:
//	C_Progress_display progress;
//
//	// Load the SfM_Data region's views
//	std::shared_ptr<Regions_Provider> regions_provider = std::make_shared<Regions_Provider>();
//	if (!regions_provider->load(sfm_data, sOutputDir, regions_type, &progress)) {
//		std::cerr << std::endl << "Invalid regions." << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	if (!stlplus::folder_exists(sQueryDir) && !stlplus::file_exists(sQueryDir))
//	{
//		std::cerr << "\nThe query directory/file does not exist : " << std::endl;
//		std::cerr << sQueryDir << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	if (sLocOutDir.empty()) {
//		std::cerr << "\nPlease provide a valid directory for the option [-o|--out_dir]." << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	if (!stlplus::folder_exists(sLocOutDir))
//		stlplus::folder_create(sLocOutDir);
//
//	if (bUseSingleIntrinsics && sfm_data.GetIntrinsics().size() != 1)
//	{
//		std::cout << "More than one intrinsics to compare to in input scene "
//			<< " => Consider intrinsics as unkown." << std::endl;
//	}
//
//	//-- Localization
//	// - init the retrieval database
//	// - Go along the sfm_data view
//	// - extract the regions of the view
//	// - try to locate the images
//	// - add the images to the sfm_data scene
//
//	std::vector<Vec3> vec_found_poses;
//
//	sfm::SfM_Localization_Single_3DTrackObservation_Database localizer;
//	if (!localizer.Init(sfm_data, *regions_provider.get()))
//	{
//		std::cerr << "Cannot initialize the SfM localizer" << std::endl;
//	}
//	// Since we have copied interesting data, release some memory
//	regions_provider.reset();
//
//	// list images from sfm_data in a vector
//	std::vector<std::string> vec_image_original(sfm_data.GetViews().size());
//	int n(-1);
//	std::generate(vec_image_original.begin(),
//		vec_image_original.end(),
//		[&n, &sfm_data]
//	{
//		n++;
//		return stlplus::filename_part(sfm_data.views.at(n)->s_Img_path);
//	});
//
//	// list images in query directory
//	std::vector<std::string> vec_image;
//
//	if (stlplus::is_file(sQueryDir))
//	{
//		vec_image.emplace_back(stlplus::filename_part(sQueryDir)); // single file
//		sQueryDir = stlplus::folder_part(sQueryDir);
//	}
//	else vec_image = stlplus::folder_files(sQueryDir); // multiple files
//
//	std::sort(vec_image.begin(), vec_image.end());
//
//	// find difference between two list of images
//	std::vector<std::string> vec_image_new;
//	std::set_difference(vec_image.cbegin(), vec_image.cend(),
//		vec_image_original.cbegin(), vec_image_original.cend(),
//		std::back_inserter(vec_image_new));
//
//	// find common root directory between images in vec_image_original and vec_images_new
//	const std::string common_root_dir = FindCommonRootDir(sfm_data.s_root_path, sQueryDir);
//
//	// check if sfm_data's root dir differs from the common root dir.
//	if (sfm_data.s_root_path != common_root_dir)
//	{
//		// in that case we have to change all the image paths from the original
//		// reconstruction
//		for (auto & view : sfm_data.GetViews())
//		{
//			view.second->s_Img_path = stlplus::create_filespec(stlplus::folder_to_relative_path(common_root_dir, sfm_data.s_root_path),
//				view.second->s_Img_path);
//		}
//		// change root path to common root path
//		sfm_data.s_root_path = common_root_dir;
//	}
//
//	// references
//	Views & views = sfm_data.views;
//	Poses & poses = sfm_data.poses;
//	Intrinsics & intrinsics = sfm_data.intrinsics;
//
//	int total_num_images = 0;
//
//#ifdef OPENMVG_USE_OPENMP
//	const unsigned int nb_max_thread = (iNumThreads == 0) ? 0 : omp_get_max_threads();
//	omp_set_num_threads(nb_max_thread);
//#pragma omp parallel for schedule(dynamic)
//#endif
//	for (int i = 0; i < static_cast<int>(vec_image_new.size()); ++i)
//	{
//		std::vector<std::string>::const_iterator iter_image = vec_image_new.begin();
//		std::advance(iter_image, i);
//
//
//		// Test if the image format is supported:
//		if (openMVG::image::GetFormat((*iter_image).c_str()) == openMVG::image::Unknown)
//		{
//			std::cerr << *iter_image << " : unknown image file format." << std::endl;
//			continue;
//		}
//
//		std::cout << "SfM::localization => try with image: " << *iter_image << std::endl;
//		std::unique_ptr<Regions> query_regions(regions_type->EmptyClone());
//		image::Image<unsigned char> imageGray;
//		{
//			const std::string sView_filename = stlplus::create_filespec(sQueryDir, *iter_image);
//			// Try to open image
//			if (!image::ReadImage(sView_filename.c_str(), &imageGray))
//			{
//				std::cerr << "Cannot open the input provided image : " << *iter_image << std::endl;
//				continue;
//			}
//
//			const std::string
//				sFeat = stlplus::create_filespec(sMatchesOutDir, stlplus::basename_part(sView_filename.c_str()), "feat"),
//				sDesc = stlplus::create_filespec(sMatchesOutDir, stlplus::basename_part(sView_filename.c_str()), "desc");
//
//			// Compute features and descriptors and save them if they don't exist yet
//			if (!stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc))
//			{
//				image_describer->Describe(imageGray, query_regions);
//				image_describer->Save(query_regions.get(), sFeat, sDesc);
//				std::cout << "#regions detected in query image: " << query_regions->RegionCount() << std::endl;
//			}
//			else // load already existing regions
//			{
//				query_regions->Load(sFeat, sDesc);
//			}
//		}
//
//		std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic;
//		if (bUseSingleIntrinsics)
//		{
//			if (sfm_data.GetIntrinsics().size() != 1)
//			{
//				std::cerr << "You choose the single intrinsic mode but the sfm_data scene,"
//					<< " have too few or too much intrinsics."
//					<< std::endl;
//				continue;
//			}
//			optional_intrinsic = sfm_data.GetIntrinsics().at(0);
//			if (imageGray.Width() != optional_intrinsic->w() || optional_intrinsic->h() != imageGray.Height())
//			{
//				std::cout << "The provided image does not have the same size as the camera model you want to use." << std::endl;
//				continue;
//			}
//		}
//		if (optional_intrinsic)
//		{
//			std::cout << "- use known intrinsics." << std::endl;
//		}
//		else
//		{
//			std::cout << "- use Unknown intrinsics for the resection. A new camera (intrinsic) will be created." << std::endl;
//
//			// Since the spherical image is only defined by its image size we can initialize its camera model.
//			// This way the resection will be performed with valid bearing vector
//			if (openMVG::cameras::EINTRINSIC(i_User_camera_model) == cameras::CAMERA_SPHERICAL)
//			{
//				optional_intrinsic = std::make_shared<cameras::Intrinsic_Spherical>(imageGray.Width(), imageGray.Height());
//			}
//		}
//
//		geometry::Pose3 pose;
//		sfm::Image_Localizer_Match_Data matching_data;
//		matching_data.error_max = dMaxResidualError;
//
//		bool bSuccessfulLocalization = false;
//
//		// Try to localize the image in the database thanks to its regions
//		if (!localizer.Localize(
//			optional_intrinsic ? resection::SolverType::P3P_KE_CVPR17 : resection::SolverType::DLT_6POINTS,
//			{ imageGray.Width(), imageGray.Height() },
//			optional_intrinsic.get(),
//			*(query_regions.get()),
//			pose,
//			&matching_data))
//		{
//			std::cerr << "Cannot locate the image " << *iter_image << std::endl;
//			bSuccessfulLocalization = false;
//		}
//		else
//		{
//			const bool b_new_intrinsic = (optional_intrinsic == nullptr);
//			// A valid pose has been found (try to refine it):
//			// If not intrinsic as input:
//			// init a new one from the projection matrix decomposition
//			// Else use the existing one and consider as static.
//			if (b_new_intrinsic)
//			{
//				// setup a default camera model from the found projection matrix
//				Mat3 K, R;
//				Vec3 t;
//				KRt_From_P(matching_data.projection_matrix, &K, &R, &t);
//
//				const double focal = (K(0, 0) + K(1, 1)) / 2.0;
//				const Vec2 principal_point(K(0, 2), K(1, 2));
//
//				switch (openMVG::cameras::EINTRINSIC(i_User_camera_model))
//				{
//				case cameras::PINHOLE_CAMERA:
//					optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic>(imageGray.Width(), imageGray.Height(), focal, principal_point(0), principal_point(1));
//					break;
//				case cameras::PINHOLE_CAMERA_RADIAL1:
//					optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Radial_K1>(imageGray.Width(), imageGray.Height(), focal, principal_point(0), principal_point(1));
//					break;
//				case cameras::PINHOLE_CAMERA_RADIAL3:
//					optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Radial_K3>(imageGray.Width(), imageGray.Height(), focal, principal_point(0), principal_point(1));
//					break;
//				case cameras::PINHOLE_CAMERA_BROWN:
//					optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Brown_T2>(imageGray.Width(), imageGray.Height(), focal, principal_point(0), principal_point(1));
//					break;
//				case cameras::PINHOLE_CAMERA_FISHEYE:
//					optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Fisheye>(imageGray.Width(), imageGray.Height(), focal, principal_point(0), principal_point(1));
//					break;
//				case cameras::CAMERA_SPHERICAL:
//					std::cerr << "The spherical camera cannot be created there. Resection of a spherical camera must be done with an existing camera model." << std::endl;
//					break;
//				default:
//					std::cerr << "Error: unknown camera model: " << static_cast<int>(i_User_camera_model) << std::endl;
//				}
//			}
//			if (optional_intrinsic && sfm::SfM_Localizer::RefinePose(
//				optional_intrinsic.get(),
//				pose, matching_data,
//				true, b_new_intrinsic))
//			{
//				bSuccessfulLocalization = true;
//			}
//			else
//			{
//				std::cerr << "Refining pose for the image " << *iter_image << " failed." << std::endl;
//			}
//
//		}
//#ifdef OPENMVG_USE_OPENMP
//#pragma omp critical
//#endif
//		{
//			total_num_images++;
//
//			View v(*iter_image, views.size(), views.size(), views.size(), imageGray.Width(), imageGray.Height());
//			if (bSuccessfulLocalization)
//			{
//				vec_found_poses.push_back(pose.center());
//				std::cout << v.s_Img_path << "," << pose.center()[0] << "," << pose.center()[1] << "," << pose.center()[2] << std::endl;
//				out<< stlplus::filename_part(v.s_Img_path)<<"," << pose.center()[0] << "," << pose.center()[1] << "," << pose.center()[2] << std::endl;
//				// Add the computed intrinsic to the sfm_container
//				if (!bUseSingleIntrinsics)
//					intrinsics[v.id_intrinsic] = optional_intrinsic;
//				else // Make the view using the existing intrinsic id
//					v.id_intrinsic = sfm_data.GetViews().begin()->second->id_intrinsic;
//				// Add the computed pose to the sfm_container
//				poses[v.id_pose] = pose;
//
//			}
//			else
//			{
//				v.id_intrinsic = UndefinedIndexT;
//				v.id_pose = UndefinedIndexT;
//			}
//			// Add the view to the sfm_container
//			views[v.id_view] = std::make_shared<View>(v);
//		}
//	}
//	
//	GroupSharedIntrinsics(sfm_data);
//	std::cout << " Total poses found : " << vec_found_poses.size() << "/" << total_num_images << endl;
//	out << ";Total poses found :" << vec_found_poses.size() << "/" << total_num_images << endl;
//	out.close();
//	return EXIT_SUCCESS;
//}

int CSFM::groupSfM()
{
	//sfm_init_MCImageListing();
	//computeMCFeatures();
	//computeMCMatches();
	globalMCSfM();
	//incrementalMCSfM2();
	Evaluate_InitialPoses();
	return EXIT_SUCCESS;
}

int main()
{
	//sfm_init_ImageListing
	CSFM sfm;
	sfm.sImageDir = "E:/qj/concentration/gopro_cc/undistorted_photos";
	sfm.sOutputDir = "D:/1graduation/suidao/";
	sfm.sSfM_Data_Filename = sfm.sOutputDir + "sfm_data.json";
	sfm.recordfile = "D:/1graduation/suidao/record.txt";
	sfm.i_User_camera_model = PINHOLE_CAMERA;
	sfm.b_use_motion_priors = false;
	sfm.Group_camera_num = 7;
	sfm.sKmatrix = "1754.5025070765;1999.15616982878;1504.93517804714;1756.32891440917;1999.69337243298;1516.10350889407;1762.1039454281;1998.12787930753;1511.8461723193;1758.94141840993;1984.91635244305;1508.77780339745;1753.47218053945;1993.51101954389;1520.61457333347;1755.63373929306;2006.60312067035;1507.26947519801;1757.83574549348;1994.89306111102;1522.69675961653;";
	sfm.ePairmode = PAIR_FOR_MULTICAMERAS;
	sfm.b_use_motion_priors = false;
	sfm.sImage_Describer_Method = "SIFT_GPU";
	sfm.MCRMSEs.resize(sfm.Group_camera_num, 0);
	sfm.groupSfM();
	return getchar();
}

