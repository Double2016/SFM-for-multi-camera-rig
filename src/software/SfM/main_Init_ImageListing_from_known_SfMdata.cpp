#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <atomic>
#include <cstdlib>

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"
#include "openMVG/exif/sensor_width_database/ParseDatabase.hpp"
#include "nonFree/sift/SIFT_describer_io.hpp"
#include "openMVG/geodesy/geodesy.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/types.hpp"

#include <cereal/archives/json.hpp>

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::cameras;
using namespace openMVG::image;
using namespace openMVG::exif;

std::pair<bool, Vec3> checkPriorWeightsString_(const std::string &sWeights)
{
	std::pair<bool, Vec3> val(true, Vec3::Zero());
	std::vector<std::string> vec_str;
	stl::split(sWeights, ';', vec_str);
	if (vec_str.size() != 3)
	{
		std::cerr << "\n Missing ';' character in prior weights" << std::endl;
		val.first = false;
	}
	// Check that all weight values are valid numbers
	for (size_t i = 0; i < vec_str.size(); ++i)
	{
		double readvalue = 0.0;
		std::stringstream ss;
		ss.str(vec_str[i]);
		if (!(ss >> readvalue)) {
			std::cerr << "\n Used an invalid not a number character in local frame origin" << std::endl;
			val.first = false;
		}
		val.second[i] = readvalue;
	}
	return val;
}
bool checkIntrinsicStringValidity_(const std::string & Kmatrix, double & focal, double & ppx, double & ppy)
{
	std::vector<std::string> vec_str;
	stl::split(Kmatrix, ';', vec_str);
	if (vec_str.size() != 9) {
		std::cerr << "\n Missing ';' character" << std::endl;
		return false;
	}
	// Check that all K matrix value are valid numbers
	for (size_t i = 0; i < vec_str.size(); ++i) {
		double readvalue = 0.0;
		std::stringstream ss;
		ss.str(vec_str[i]);
		if (!(ss >> readvalue)) {
			std::cerr << "\n Used an invalid not a number character" << std::endl;
			return false;
		}
		if (i == 0) focal = readvalue;
		if (i == 2) ppx = readvalue;
		if (i == 5) ppy = readvalue;
	}
	return true;
}
int main()
{
	SfM_Data m_sfm_data;
	int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;
	std::string sImageDir = "E:/panxian/photo";
	std::string sOutputDir = "E:/panxian/test_ppose";
	std::string SensorDatafile = "D:/openMVG-develop/openMVG-develop/src/openMVG/exif/sensor_width_database/sensor_width_camera_database.txt";
	std::string sfm_data_BA_File = "E:/panxian/reconstruction_global/sfm_data_vie.json";
	bool b_Group_camera_model = true;
	double focal_pixels = -1.0;
	std::string sKmatrix = "";
	bool use_pose_prior = true;
	std::string sPriorWeights;

	//CmdLine cmd;
	std::cout << "--------------------------sfm_init_ImageListing-----------------------------" << std::endl;
	std::pair<bool, Vec3> prior_w_info(false, Vec3(1.0, 1.0, 1.0));

	std::cout << " You called : " << std::endl
		<< "--imageDirectory " << sImageDir << std::endl
		<< "--sensorWidthDatabase " << SensorDatafile << std::endl
		<< "--outputDirectory " << sOutputDir << std::endl
		<< "--focal " << focal_pixels << std::endl
		<< "--intrinsics " << sKmatrix << std::endl
		<< "--camera_model " << i_User_camera_model << std::endl
		<< "--group_camera_model " << b_Group_camera_model << std::endl
		<< "-- use_pose_prior " << use_pose_prior << std::endl;

	// Expected properties for each image
	double width = -1, height = -1, focal = -1, ppx = -1, ppy = -1;

	const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

	if (!stlplus::folder_exists(sImageDir))
	{
		std::cerr << "\nThe input directory doesn't exist" << std::endl;
		return getchar();
	}

	if (sOutputDir.empty())
	{
		std::cerr << "\nInvalid output directory" << std::endl;
		return getchar();
	}

	if (!stlplus::folder_exists(sOutputDir))
	{
		if (!stlplus::folder_create(sOutputDir))
		{
			std::cerr << "\nCannot create output directory" << std::endl;
			return getchar();
		}
	}

	if (sKmatrix.size() > 0 &&
		!checkIntrinsicStringValidity_(sKmatrix, focal, ppx, ppy))
	{
		std::cerr << "\nInvalid K matrix input" << std::endl;
		return getchar();
	}

	if (sKmatrix.size() > 0 && focal_pixels != -1.0)
	{
		std::cerr << "\nCannot combine -f and -k options" << std::endl;
		return getchar();
	}

	std::vector<Datasheet> vec_database;
	if (!SensorDatafile.empty())
	{
		if (!parseDatabase(SensorDatafile, vec_database))
		{
			std::cerr
				<< "\nInvalid input database: " << SensorDatafile
				<< ", please specify a valid file." << std::endl;
			return getchar();
		}
	}

	// Check if prior weights are given
	if (use_pose_prior && !sPriorWeights.empty())
	{
		prior_w_info = checkPriorWeightsString_(sPriorWeights);
	}
	else if (use_pose_prior)
	{
		prior_w_info.first = true;
	}

	std::vector<std::string> vec_image = stlplus::folder_files(sImageDir);
	std::sort(vec_image.begin(), vec_image.end());

	//Verified by Double
	SfM_Data ref_sfm_data;
	Load(ref_sfm_data, sfm_data_BA_File, ESfM_Data(VIEWS | EXTRINSICS|INTRINSICS));

	if (ref_sfm_data.views.size() == vec_image.size())
	{
		std::cout << "Begin to use a referenced sfm_data as prior pose...." << std::endl;
		std::cout << "ref_sfm_data.poses.size()=" << ref_sfm_data.poses.size() << std::endl;
	}
	else
	{
		std::cerr << "\nFailed!!!" << std::endl;
		return getchar();
	}

	//Verified by Double

	// Configure an empty scene with Views and their corresponding cameras
	SfM_Data sfm_data;
	sfm_data.s_root_path = sImageDir; // Setup main image root_path
	Views & views = sfm_data.views;
	Intrinsics & intrinsics = sfm_data.intrinsics;

	//C_Progress_display my_progress_bar(vec_image.size(),
	//	std::cout, "\n- Image listing -\n");
	//std::ostringstream error_report_stream;
	for (const auto & view_it : ref_sfm_data.GetViews())
	{
		// Read meta data to fill camera parameter (w,h,focal,ppx,ppy) fields.
		width = height = ppx = ppy = focal = -1.0;

		const std::string sImageFilename = stlplus::create_filespec(sImageDir,view_it.second->s_Img_path);
		const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

		// Test if the image format is supported:
		if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown)
		{
			std::cerr
				<< sImFilenamePart << ": Unkown image file format." << "\n";
			continue; // image cannot be opened
		}

		if (sImFilenamePart.find("mask.png") != std::string::npos
			|| sImFilenamePart.find("_mask.png") != std::string::npos)
		{
			std::cerr
				<< sImFilenamePart << " is a mask image" << "\n";
			continue;
		}

		image::ImageHeader imgHeader;
		if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
			continue; // image cannot be read

		width = imgHeader.width;
		height = imgHeader.height;
		ppx = width / 2.0;
		ppy = height / 2.0;


		// Consider the case where the focal is provided manually
		if (sKmatrix.size() > 0) // Known user calibration K matrix
		{
			if (!checkIntrinsicStringValidity_(sKmatrix, focal, ppx, ppy))
				focal = -1.0;
		}
		else // User provided focal length value
			if (focal_pixels != -1)
				focal = focal_pixels;

		// If not manually provided or wrongly provided
		if (focal == -1)
		{
			std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
			exifReader->open(sImageFilename);

			const bool bHaveValidExifMetadata =
				exifReader->doesHaveExifInfo()
				&& !exifReader->getModel().empty();

			if (bHaveValidExifMetadata) // If image contains meta data
			{
				const std::string sCamModel = exifReader->getModel();

				// Handle case where focal length is equal to 0
				if (exifReader->getFocal() == 0.0f)
				{
					std::cerr
						<< stlplus::basename_part(sImageFilename) << ": Focal length is missing." << "\n";
					focal = -1.0;
				}
				else
					// Create the image entry in the list file
				{
					Datasheet datasheet;
					if (getInfo(sCamModel, vec_database, datasheet))
					{
						// The camera model was found in the database so we can compute it's approximated focal length
						const double ccdw = datasheet.sensorSize_;
						focal = (std::max)(width, height) * exifReader->getFocal() / ccdw;
					}
					else
					{
						std::cerr
							<< stlplus::basename_part(sImageFilename)
							<< "\" model \"" << sCamModel << "\" doesn't exist in the database" << "\n"
							<< "Please consider add your camera model and sensor width in the database." << "\n";
					}
				}
			}
		}
		// Build intrinsic parameter related to the view
		std::shared_ptr<IntrinsicBase> intrinsic;

		if (focal > 0 && ppx > 0 && ppy > 0 && width > 0 && height > 0)
		{
			// Create the desired camera type
			switch (e_User_camera_model)
			{
			case PINHOLE_CAMERA:
				intrinsic = std::make_shared<Pinhole_Intrinsic>
					(width, height, focal, ppx, ppy);
				break;
			case PINHOLE_CAMERA_RADIAL1:
				intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
					(width, height, focal, ppx, ppy, 0.0); // setup no distortion as initial guess
				break;
			case PINHOLE_CAMERA_RADIAL3:
				intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
					(width, height, focal, ppx, ppy, 0.0, 0.0, 0.0);  // setup no distortion as initial guess
				break;
			case PINHOLE_CAMERA_BROWN:
				intrinsic = std::make_shared<Pinhole_Intrinsic_Brown_T2>
					(width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
				break;
			case PINHOLE_CAMERA_FISHEYE:
				intrinsic = std::make_shared<Pinhole_Intrinsic_Fisheye>
					(width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
				break;
			case CAMERA_SPHERICAL:
				intrinsic = std::make_shared<Intrinsic_Spherical>
					(width, height);
				break;
			default:
				std::cerr << "Error: unknown camera model: " << (int)e_User_camera_model << std::endl;
				return getchar();
			}
		}

		// Build the view corresponding to the image

		std::pair<bool, Vec3> gps_info;//= checkGPSfromfile(path_PoseFile, sImFilenamePart);

		if ( use_pose_prior)
		{
			ViewPriors v(view_it.second->s_Img_path, views.size(), views.size(), views.size(), width, height);
			if (ref_sfm_data.IsPoseAndIntrinsicDefined(view_it.second.get()))
			{
				gps_info.first = true;
				gps_info.second = ref_sfm_data.GetPoses().at(view_it.second->id_pose).center();
			}
			else
			{
				std::cout << "Due to this view's pose is not Defined,Failed!!!!!" << std::endl;
				continue;
			}
				
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
			if (gps_info.first)
			{
				v.b_use_pose_center_ = true;
				v.pose_center_ = gps_info.second;
				// prior weights
				if (prior_w_info.first == true)
				{
					v.center_weight_ = prior_w_info.second;
				}
				// Add the view to the sfm_container
				views[v.id_view] = std::make_shared<ViewPriors>(v);
			}
			
		}
		else
		{
			std::cerr << "\n Failed!........." << std::endl;
			return getchar();
		}
	}


	// Group camera that share common properties if desired (leads to more faster & stable BA).
	if (b_Group_camera_model)
	{
		GroupSharedIntrinsics(sfm_data);
	}

	if (!Save(
		sfm_data,
		stlplus::create_filespec(sOutputDir, "sfm_data.json").c_str(),
		ESfM_Data(VIEWS | INTRINSICS | CONTROL_POINTS))) // VIEWS | INTRINSIC ΪɶBUG
	{
		return getchar();
	}

	std::cout << std::endl
		<< "SfMInit_ImageListing report:\n"
		<< "listed #File(s): " << vec_image.size() << "\n"
		<< "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
		<< "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size() << std::endl;

	return getchar();

}