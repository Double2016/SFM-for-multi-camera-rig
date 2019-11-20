#ifndef SIFT_GPU_H
#define SIFT_GPU_H

#include<iostream>

void hello_siftgpu();


#include<iostream>
#include <cereal/archives/json.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/polymorphic.hpp>
#include "openMVG/image/image_io.hpp"
#include "openMVG/features/regions_factory_io.hpp"
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/system/timer.hpp"


//#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"
#define USE_OCVSIFTs
#define USE_SIFTGPU
/// OpenCV Includes
//#include <opencv2/opencv.hpp>
//#include "opencv2/core/eigen.hpp"
#ifdef USE_OCVSIFT
#include "opencv2/xfeatures2d.hpp"
#endif

#include "SiftGPU.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::features;



//namespace openMVG {
//	namespace features {
//		//qu modify add start
//		using SIFT_GPU_Regions = Scalar_Regions<SIOPointFeature, float, 128>;
//		//qu modify add end
//	} // namespace features
//} // namespace openMVG
//
//  //qu modify add start
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST(openMVG::features::SIFT_GPU_Regions)
////qu modify add end
//
//
////q modify add start
//CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::SIFT_GPU_Regions, "SIFT_GPU_Regions");
//CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Regions, openMVG::features::SIFT_GPU_Regions)
////q modify add end




#ifdef USE_OCVSIFT
///
//- Create an Image_describer interface that use and OpenCV extraction method
// i.e. with the SIFT detector+descriptor
// Regions is the same as classic SIFT : 128 unsigned char
class SIFT_OPENCV_Image_describer : public Image_describer
{
public:
	using Regions_type = SIFT_Regions;

	SIFT_OPENCV_Image_describer() : Image_describer() {}

	~SIFT_OPENCV_Image_describer() {}

	bool Set_configuration_preset(EDESCRIBER_PRESET preset) {
		return true;
	}

	/**
	@brief Detect regions on the image and compute their attributes (description)
	@param image Image.
	@param mask 8-bit gray image for keypoint filtering (optional).
	Non-zero values depict the region of interest.
	@return regions The detected regions and attributes (the caller must delete the allocated data)
	*/
	std::unique_ptr<Regions> Describe(
		const image::Image<unsigned char>& image,
		const image::Image<unsigned char> * mask = nullptr
	) override
	{
		return Describe_SIFT_OPENCV(image, mask);
	}

	/**
	@brief Detect regions on the image and compute their attributes (description)
	@param image Image.
	@param mask 8-bit gray image for keypoint filtering (optional).
	Non-zero values depict the region of interest.
	@return regions The detected regions and attributes (the caller must delete the allocated data)
	*/
	std::unique_ptr<Regions_type> Describe_SIFT_OPENCV(
		const image::Image<unsigned char>& image,
		const image::Image<unsigned char>* mask = nullptr
	)
	{
		// Convert for opencv
		cv::Mat img;
		cv::eigen2cv(image.GetMat(), img);

		// Convert mask image into cv::Mat
		cv::Mat m_mask;
		if (mask != nullptr) {
			cv::eigen2cv(mask->GetMat(), m_mask);
		}

		// Create a SIFT detector
		std::vector< cv::KeyPoint > v_keypoints;
		cv::Mat m_desc;
		cv::Ptr<cv::Feature2D> siftdetector = cv::xfeatures2d::SIFT::create();

		// Process SIFT computation
		siftdetector->detectAndCompute(img, m_mask, v_keypoints, m_desc);

		auto regions = std::unique_ptr<Regions_type>(new Regions_type);

		// reserve some memory for faster keypoint saving
		regions->Features().reserve(v_keypoints.size());
		regions->Descriptors().reserve(v_keypoints.size());

		// Prepare a column vector with the sum of each descriptor
		cv::Mat m_siftsum;
		cv::reduce(m_desc, m_siftsum, 1, cv::REDUCE_SUM);

		// Copy keypoints and descriptors in the regions
		int cpt = 0;
		for (auto i_kp = v_keypoints.begin();
			i_kp != v_keypoints.end();
			++i_kp, ++cpt)
		{
			SIOPointFeature feat((*i_kp).pt.x, (*i_kp).pt.y, (*i_kp).size, (*i_kp).angle);
			regions->Features().push_back(feat);

			Descriptor<unsigned char, 128> desc;
			for (int j = 0; j < 128; j++)
			{
				desc[j] = static_cast<unsigned char>(512.0*sqrt(m_desc.at<float>(cpt, j) / m_siftsum.at<float>(cpt, 0)));
			}
			regions->Descriptors().push_back(desc);
		}

		return regions;
	};

	/// Allocate Regions type depending of the Image_describer
	std::unique_ptr<Regions> Allocate() const override
	{
		return std::unique_ptr<Regions_type>(new Regions_type);
	}

	template<class Archive>
	void serialize(Archive & ar)
	{
	}
};
CEREAL_REGISTER_TYPE_WITH_NAME(SIFT_OPENCV_Image_describer, "SIFT_OPENCV_Image_describer");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, SIFT_OPENCV_Image_describer)
#endif //USE_OCVSIFT



#ifdef USE_SIFTGPU
///
/******************* use opencv ********************/
//- Create an Image_describer interface that use and OpenCV extraction method
// i.e. with the SIFT detector+descriptor
// Regions is the same as classic SIFT : 128 unsigned char
//class SIFT_GPU_Image_describer : public Image_describer
//{
//public:
//	using Regions_type = SIFT_GPU_Regions;
//
//	SIFT_GPU_Image_describer() : Image_describer() 
//	{
//		char * argv1[] = { "-fo", "0",  "-v", "0","-maxd","19200", "-t","0.016" };
//		int argc1 = sizeof(argv1) / sizeof(char*);
//
//		sift.ParseParam(argc1, argv1);
//		if (sift.CreateContextGL() != SiftGPU::SIFTGPU_FULL_SUPPORTED)
//		{
//			std::cout << "siftgpu not support !" << std::endl;
//			return;
//		}
//		//sift.AllocatePyramid(2832, 2128);
//	}
//
//	~SIFT_GPU_Image_describer() {}
//
//	bool Set_configuration_preset(EDESCRIBER_PRESET preset) {
//		return true;
//	}
//
//	void initSiftGPU(int _width, int _height)
//	{
//		width = _width;
//		height = _height;
//		//Maximum working dimension. When some level images are larger
//		//than this, the input image will be automatically down - sampled.
//		//int maxTextureDimension = max(_width, _height);
//		sift.AllocatePyramid(_width, _height);
//	}
//
//	/**
//	@brief Detect regions on the image and compute their attributes (description)
//	@param image Image.
//	@param mask 8-bit gray image for keypoint filtering (optional).
//	Non-zero values depict the region of interest.
//	@return regions The detected regions and attributes (the caller must delete the allocated data)
//	*/
//	std::unique_ptr<Regions> Describe(
//		const image::Image<unsigned char>& image,
//		const image::Image<unsigned char> * mask = nullptr
//	) override
//	{
//		return Describe_SIFT_OPENCV(image, mask);
//	}
//
//	/**
//	@brief Detect regions on the image and compute their attributes (description)
//	@param image Image.
//	@param mask 8-bit gray image for keypoint filtering (optional).
//	Non-zero values depict the region of interest.
//	@return regions The detected regions and attributes (the caller must delete the allocated data)
//	*/
//	std::unique_ptr<Regions_type> Describe_SIFT_OPENCV(
//		const image::Image<unsigned char>& image,
//		const image::Image<unsigned char>* mask = nullptr
//	)
//	{
//		// Convert for opencv
//		cv::Mat img;
//		cv::eigen2cv(image.GetMat(), img);
//
//		// Convert mask image into cv::Mat
//		cv::Mat m_mask;
//		if (mask != nullptr) {
//			cv::eigen2cv(mask->GetMat(), m_mask);
//		}
//
//		//// Create a SIFT detector
//		//std::vector< cv::KeyPoint > v_keypoints;
//		//cv::Mat m_desc;
//		//cv::Ptr<cv::Feature2D> siftdetector = cv::xfeatures2d::SIFT::create();
//
//		//// Process SIFT computation
//		//siftdetector->detectAndCompute(img, m_mask, v_keypoints, m_desc);
//		
//		//char * argv1[] = { "-fo", "0",  "-v", "0" };//
//
//		using namespace cv;
//		cv::Mat m_desc;
//		std::vector<SiftGPU::SiftKeypoint> v_keypoints;
//		int width = img.cols;
//		int height = img.rows;
//		int num;
//		//GL_LUMINANCE 0x1909 GL_UNSIGNED_BYTE 0x1401
//		// (gl_format == GL_LUMINANCE || gl_format == GL_LUMINANCE_ALPHA ||
//		//	gl_format == GL_RGB || gl_format == GL_RGBA ||
//		//	gl_format == GL_BGR || gl_format == GL_BGRA) &&
//		//	(gl_type == GL_UNSIGNED_BYTE || gl_type == GL_FLOAT || gl_type == GL_UNSIGNED_SHORT);
//		if (sift.RunSIFT(width, height, img.data, 0x1909, 0x1401))
//		{
//			//Call SaveSIFT to save result to file, the format is the same as Lowe's
//			//sift->SaveSIFT("../data/800-1.sift"); //Note that saving ASCII format is slow
//
//			//get feature count
//			num = sift.GetFeatureNum();
//
//			//allocate memory
//			v_keypoints.resize(num);    m_desc.create(num,128,CV_32F);
//
//			//reading back feature vectors is faster than writing files
//			//if you dont need keys or descriptors, just put NULLs here
//			sift.GetFeatureVector(&v_keypoints[0], (float*)m_desc.data);
//			//this can be used to write your own sift file.            
//		}
//
//		auto regions = std::unique_ptr<Regions_type>(new Regions_type);
//
//		// reserve some memory for faster keypoint saving
//		regions->Features().reserve(v_keypoints.size());
//		regions->Descriptors().reserve(v_keypoints.size());
//
//		// Prepare a column vector with the sum of each descriptor
//		//cv::Mat m_siftsum;
//		//cv::reduce(m_desc, m_siftsum, 1, cv::REDUCE_SUM);
//
//		// Copy keypoints and descriptors in the regions
//		//int cpt = 0;
//		//for (auto i_kp = v_keypoints.begin();
//		//	i_kp != v_keypoints.end();
//		//	++i_kp, ++cpt)
//		//{
//		//	SIOPointFeature feat((*i_kp).x, (*i_kp).y, (*i_kp).s, (*i_kp).o);
//		//	regions->Features().push_back(feat);
//
//		//	Descriptor<unsigned char, 128> desc;
//		//	const float* p_des = m_desc.ptr<float>(cpt);
//		//	for (int j = 0; j < 128; j++)
//		//	{
//		//		desc[j] = static_cast<unsigned char>(512.0*sqrt(m_desc.at<float>(cpt, j) / m_siftsum.at<float>(cpt, 0)));
//		//	}
//		//	regions->Descriptors().push_back(desc);
//		//}
//		int cpt = 0;
//		for (auto i_kp = v_keypoints.begin();
//			i_kp != v_keypoints.end();
//			++i_kp, ++cpt)
//		{
//			SIOPointFeature feat((*i_kp).x, (*i_kp).y, (*i_kp).s, (*i_kp).o);
//			regions->Features().push_back(feat);
//
//			Descriptor<float, 128> desc;
//			const float* p_des = m_desc.ptr<float>(cpt);
//			//const float* p_sum = m_siftsum.ptr<float>(cpt);
//			for (int j = 0; j < 128; j++)
//			{
//				//desc[j] = static_cast<unsigned char>(512.0*sqrt(p_des[j] / p_sum[0]));
//				desc[j] = p_des[j];
//			}
//			regions->Descriptors().push_back(desc);
//		}
//		return regions;
//	};
//
//	/// Allocate Regions type depending of the Image_describer
//	std::unique_ptr<Regions> Allocate() const override
//	{
//		return std::unique_ptr<Regions_type>(new Regions_type);
//	}
//
//	template<class Archive>
//	void serialize(Archive & ar)
//	{
//	}
//private:
//	SiftGPU sift;
//	int width, height;
//};
/******************* no use opencv ********************/
//class SIFT_GPU_Image_describer : public Image_describer
//{
//public:
//	using Regions_type = SIFT_GPU_Regions;
//
//	SIFT_GPU_Image_describer() : Image_describer()
//	{
//		char * argv1[] = { "-fo", "0",  "-v", "0","-maxd","19200", "-t","0.016" };
//		int argc1 = sizeof(argv1) / sizeof(char*);
//
//		sift.ParseParam(argc1, argv1);
//		if (sift.CreateContextGL() != SiftGPU::SIFTGPU_FULL_SUPPORTED)
//		{
//			std::cout << "siftgpu not support !" << std::endl;
//			return;
//		}
//		//sift.AllocatePyramid(2832, 2128);
//	}
//
//	~SIFT_GPU_Image_describer() {}
//
//	bool Set_configuration_preset(EDESCRIBER_PRESET preset) {
//		return true;
//	}
//
//	void initSiftGPU(int _width, int _height)
//	{
//		width = _width;
//		height = _height;
//		//Maximum working dimension. When some level images are larger
//		//than this, the input image will be automatically down - sampled.
//		//int maxTextureDimension = max(_width, _height);
//		sift.AllocatePyramid(_width, _height);
//	}
//
//	/**
//	@brief Detect regions on the image and compute their attributes (description)
//	@param image Image.
//	@param mask 8-bit gray image for keypoint filtering (optional).
//	Non-zero values depict the region of interest.
//	@return regions The detected regions and attributes (the caller must delete the allocated data)
//	*/
//	std::unique_ptr<Regions> Describe(
//		const image::Image<unsigned char>& image,
//		const image::Image<unsigned char> * mask = nullptr
//	) override
//	{
//		return Describe_SIFT_OPENCV(image, mask);
//	}
//
//	/**
//	@brief Detect regions on the image and compute their attributes (description)
//	@param image Image.
//	@param mask 8-bit gray image for keypoint filtering (optional).
//	Non-zero values depict the region of interest.
//	@return regions The detected regions and attributes (the caller must delete the allocated data)
//	*/
//	std::unique_ptr<Regions_type> Describe_SIFT_OPENCV(
//		const image::Image<unsigned char>& image,
//		const image::Image<unsigned char>* mask = nullptr
//	)
//	{
//		// Convert for opencv
//		//cv::Mat img;
//		//cv::eigen2cv(image.GetMat(), img);
//		auto* data = image.GetMat().data();
//		// Convert mask image into cv::Mat
//		//cv::Mat m_mask;
//		//if (mask != nullptr) {
//		//	cv::eigen2cv(mask->GetMat(), m_mask);
//		//}
//
//		//// Create a SIFT detector
//		//std::vector< cv::KeyPoint > v_keypoints;
//		//cv::Mat m_desc;
//		//cv::Ptr<cv::Feature2D> siftdetector = cv::xfeatures2d::SIFT::create();
//
//		//// Process SIFT computation
//		//siftdetector->detectAndCompute(img, m_mask, v_keypoints, m_desc);
//
//		//char * argv1[] = { "-fo", "0",  "-v", "0" };//
//
//		//using namespace cv;
//		//cv::Mat m_desc;
//		float* m_desc;
//		std::vector<SiftGPU::SiftKeypoint> v_keypoints;
//		int width = image.Width();
//		int height = image.Height();
//		int num;
//		//GL_LUMINANCE 0x1909 GL_UNSIGNED_BYTE 0x1401
//		// (gl_format == GL_LUMINANCE || gl_format == GL_LUMINANCE_ALPHA ||
//		//	gl_format == GL_RGB || gl_format == GL_RGBA ||
//		//	gl_format == GL_BGR || gl_format == GL_BGRA) &&
//		//	(gl_type == GL_UNSIGNED_BYTE || gl_type == GL_FLOAT || gl_type == GL_UNSIGNED_SHORT);
//		if (sift.RunSIFT(width, height, data, 0x1909, 0x1401))
//		{
//			//Call SaveSIFT to save result to file, the format is the same as Lowe's
//			//sift->SaveSIFT("../data/800-1.sift"); //Note that saving ASCII format is slow
//
//			//get feature count
//			num = sift.GetFeatureNum();
//
//			//allocate memory
//			v_keypoints.resize(num);    
//			//m_desc.create(num, 128, CV_32F);
//			m_desc = new float[num * 128];
//
//			//reading back feature vectors is faster than writing files
//			//if you dont need keys or descriptors, just put NULLs here
//			sift.GetFeatureVector(&v_keypoints[0], (float*)m_desc);
//			//this can be used to write your own sift file.            
//		}
//
//		auto regions = std::unique_ptr<Regions_type>(new Regions_type);
//
//		// reserve some memory for faster keypoint saving
//		regions->Features().reserve(v_keypoints.size());
//		regions->Descriptors().reserve(v_keypoints.size());
//
//		// Prepare a column vector with the sum of each descriptor
//		//cv::Mat m_siftsum;
//		//cv::reduce(m_desc, m_siftsum, 1, cv::REDUCE_SUM);
//
//		// Copy keypoints and descriptors in the regions
//		//int cpt = 0;
//		//for (auto i_kp = v_keypoints.begin();
//		//	i_kp != v_keypoints.end();
//		//	++i_kp, ++cpt)
//		//{
//		//	SIOPointFeature feat((*i_kp).x, (*i_kp).y, (*i_kp).s, (*i_kp).o);
//		//	regions->Features().push_back(feat);
//
//		//	Descriptor<unsigned char, 128> desc;
//		//	const float* p_des = m_desc.ptr<float>(cpt);
//		//	for (int j = 0; j < 128; j++)
//		//	{
//		//		desc[j] = static_cast<unsigned char>(512.0*sqrt(m_desc.at<float>(cpt, j) / m_siftsum.at<float>(cpt, 0)));
//		//	}
//		//	regions->Descriptors().push_back(desc);
//		//}
//		int cpt = 0;
//		for (auto i_kp = v_keypoints.begin();
//			i_kp != v_keypoints.end();
//			++i_kp, ++cpt)
//		{
//			SIOPointFeature feat((*i_kp).x, (*i_kp).y, (*i_kp).s, (*i_kp).o);
//			regions->Features().push_back(feat);
//
//			Descriptor<float, 128> desc;
//			//const float* p_des = m_desc.ptr<float>(cpt);
//			//const float* p_sum = m_siftsum.ptr<float>(cpt);
//
//			for (int j = 0; j < 128; j++)
//			{
//				//desc[j] = static_cast<unsigned char>(512.0*sqrt(p_des[j] / p_sum[0]));
//				desc[j] = m_desc[cpt*128 + j];
//			}
//			regions->Descriptors().push_back(desc);
//		}
//		return regions;
//	};
//
//	/// Allocate Regions type depending of the Image_describer
//	std::unique_ptr<Regions> Allocate() const override
//	{
//		return std::unique_ptr<Regions_type>(new Regions_type);
//	}
//
//	template<class Archive>
//	void serialize(Archive & ar)
//	{
//	}
//private:
//	SiftGPU sift;
//	int width, height;
//};
//CEREAL_REGISTER_TYPE_WITH_NAME(SIFT_GPU_Image_describer, "SIFT_GPU_Image_describer");
//CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, SIFT_GPU_Image_describer)

/******************* no use opencv uchar descriptor ********************/
class SIFT_GPU_Image_describer : public Image_describer
{
public:
	using Regions_type = SIFT_Regions;
	struct Params
	{
		Params(
			//int first_octave = 0,
			//int num_octaves = 6,
			//int num_scales = 3,
			//float edge_threshold = 10.0f,
			double t = 0.0125
		) :
			//_first_octave(first_octave),
			//_num_octaves(num_octaves),
			//_num_scales(num_scales),
			//_edge_threshold(edge_threshold),
			_t(t)
		{}

		template<class Archive>
		inline void serialize(Archive & ar);

		// Parameters
		//int _first_octave;      // Use original image, or perform an upscale if == -1
		//int _num_octaves;       // Max octaves count
		//int _num_scales;        // Scales per octave
		//float _edge_threshold;  // Max ratio of Hessian eigenvalues
		double _t;              // DOG threshold (default : 0.0125)
	};

	SIFT_GPU_Image_describer(
		const Params & params = Params()
	) : Image_describer(), params_(params)
	{
		double t = params._t;
		//t的值不能太大，也不能太小0.02到0.02/3之间好，t越大，特征点越少（理论上精选特征点）
		if (t < 0.02 / 3)
			t = 0.02 / 3;
		if (t > 0.02)
			t = 0.02;
		//	<<"-t <float>        : DOG threshold (default : 0.02/3)\n"
		//	<<"-maxd <int> *     : Max working dimension (default : 2560 (unpacked) / 3200 (packed))\n"

		std::stringstream stream;
		stream << t;
		std::string t_str = stream.str();
		const char* c_char = t_str.c_str();
		char* t_char = const_cast<char*>(c_char);

		//char * argv1[] = { "-fo", "0",  "-v", "0","-maxd","19200", "-t","0.016" };
		char * argv1[] = { "-fo", "0",  "-v", "0","-maxd","19200", "-t",t_char };
		int argc1 = sizeof(argv1) / sizeof(char*);

		sift.ParseParam(argc1, argv1);
		if (sift.CreateContextGL() != SiftGPU::SIFTGPU_FULL_SUPPORTED)
		{
			std::cout << "siftgpu not support !" << std::endl;
			return;
		}
		//sift.AllocatePyramid(2832, 2128);
	}

	~SIFT_GPU_Image_describer() {}




	bool Set_configuration_preset(EDESCRIBER_PRESET preset) {
		return true;
	}

	void initSiftGPU(int _width, int _height)
	{
		width = _width;
		height = _height;
		//Maximum working dimension. When some level images are larger
		//than this, the input image will be automatically down - sampled.
		//int maxTextureDimension = max(_width, _height);
		sift.AllocatePyramid(_width, _height);
	}

	/**
	@brief Detect regions on the image and compute their attributes (description)
	@param image Image.
	@param mask 8-bit gray image for keypoint filtering (optional).
	Non-zero values depict the region of interest.
	@return regions The detected regions and attributes (the caller must delete the allocated data)
	*/
	std::unique_ptr<Regions> Describe(
		const image::Image<unsigned char>& image,
		const image::Image<unsigned char> * mask = nullptr
	) override
	{
		return Describe_SIFT_GPU(image, mask);
	}

	/**
	@brief Detect regions on the image and compute their attributes (description)
	@param image Image.
	@param mask 8-bit gray image for keypoint filtering (optional).
	Non-zero values depict the region of interest.
	@return regions The detected regions and attributes (the caller must delete the allocated data)
	*/
	std::unique_ptr<Regions_type> Describe_SIFT_GPU(
		const image::Image<unsigned char>& image,
		const image::Image<unsigned char>* mask = nullptr
	)
	{
		// Convert for opencv
		//cv::Mat img;
		//cv::eigen2cv(image.GetMat(), img);
		auto* data = image.GetMat().data();
		// Convert mask image into cv::Mat
		//cv::Mat m_mask;
		//if (mask != nullptr) {
		//	cv::eigen2cv(mask->GetMat(), m_mask);
		//}

		//// Create a SIFT detector
		//std::vector< cv::KeyPoint > v_keypoints;
		//cv::Mat m_desc;
		//cv::Ptr<cv::Feature2D> siftdetector = cv::xfeatures2d::SIFT::create();

		//// Process SIFT computation
		//siftdetector->detectAndCompute(img, m_mask, v_keypoints, m_desc);

		//char * argv1[] = { "-fo", "0",  "-v", "0" };//

		//using namespace cv;
		//cv::Mat m_desc;
		//float* m_desc;
		//使用智能指针不用担心其释放问题
		std::shared_ptr<float> m_desc;
		std::vector<SiftGPU::SiftKeypoint> v_keypoints;
		int width = image.Width();
		int height = image.Height();
		int num;
		//GL_LUMINANCE 0x1909 GL_UNSIGNED_BYTE 0x1401
		// (gl_format == GL_LUMINANCE || gl_format == GL_LUMINANCE_ALPHA ||
		//	gl_format == GL_RGB || gl_format == GL_RGBA ||
		//	gl_format == GL_BGR || gl_format == GL_BGRA) &&
		//	(gl_type == GL_UNSIGNED_BYTE || gl_type == GL_FLOAT || gl_type == GL_UNSIGNED_SHORT);
		if (sift.RunSIFT(width, height, data, 0x1909, 0x1401))
		{
			//Call SaveSIFT to save result to file, the format is the same as Lowe's
			//sift->SaveSIFT("../data/800-1.sift"); //Note that saving ASCII format is slow

			//get feature count
			num = sift.GetFeatureNum();

			//allocate memory
			v_keypoints.resize(num);
			//m_desc.create(num, 128, CV_32F);
			m_desc = std::shared_ptr<float>(new float[num * 128]);

			//reading back feature vectors is faster than writing files
			//if you dont need keys or descriptors, just put NULLs here
			sift.GetFeatureVector(&v_keypoints[0], (float*)m_desc.get());
			//this can be used to write your own sift file.            
		}

		auto regions = std::unique_ptr<Regions_type>(new Regions_type);

		// reserve some memory for faster keypoint saving
		regions->Features().reserve(v_keypoints.size());
		regions->Descriptors().reserve(v_keypoints.size());

		std::shared_ptr<float>m_siftsum(new float[num]);
		for (int i = 0; i < num; ++i)
		{
			float sum = .0;
			for (int j = 0; j < 128; ++j)
			{
				sum += m_desc.get()[i * 128 + j];
			}
			m_siftsum.get()[i] = sum;
		}
		// Prepare a column vector with the sum of each descriptor
		//cv::Mat m_siftsum;
		//cv::reduce(m_desc, m_siftsum, 1, cv::REDUCE_SUM);

		// Copy keypoints and descriptors in the regions
		//int cpt = 0;
		//for (auto i_kp = v_keypoints.begin();
		//	i_kp != v_keypoints.end();
		//	++i_kp, ++cpt)
		//{
		//	SIOPointFeature feat((*i_kp).x, (*i_kp).y, (*i_kp).s, (*i_kp).o);
		//	regions->Features().push_back(feat);

		//	Descriptor<unsigned char, 128> desc;
		//	const float* p_des = m_desc.ptr<float>(cpt);
		//	for (int j = 0; j < 128; j++)
		//	{
		//		desc[j] = static_cast<unsigned char>(512.0*sqrt(m_desc.at<float>(cpt, j) / m_siftsum.at<float>(cpt, 0)));
		//	}
		//	regions->Descriptors().push_back(desc);
		//}
		int cpt = 0;
		for (auto i_kp = v_keypoints.begin();
			i_kp != v_keypoints.end();
			++i_kp, ++cpt)
		{
			SIOPointFeature feat((*i_kp).x, (*i_kp).y, (*i_kp).s, (*i_kp).o);
			regions->Features().push_back(feat);

			Descriptor<unsigned char, 128> desc;
			//const float* p_des = m_desc.ptr<float>(cpt);
			//const float* p_sum = m_siftsum.ptr<float>(cpt);
			for (int j = 0; j < 128; ++j)
			{
				desc[j] = static_cast<unsigned char>(512.0*sqrt(m_desc.get()[cpt * 128 + j] / m_siftsum.get()[cpt]));
			}
			regions->Descriptors().push_back(desc);
		}
		return regions;
	};

	/// Allocate Regions type depending of the Image_describer
	std::unique_ptr<Regions> Allocate() const override
	{
		return std::unique_ptr<Regions_type>(new Regions_type);
	}

	template<class Archive>
	inline void serialize(Archive & ar);

private:
	SiftGPU sift;
	int width, height;
	Params params_;
};

template<class Archive>
inline void SIFT_GPU_Image_describer::Params::serialize(Archive & ar)
{
	ar(
		//cereal::make_nvp("first_octave", _first_octave),
		//cereal::make_nvp("num_octaves", _num_octaves),
		//cereal::make_nvp("num_scales", _num_scales),
		//cereal::make_nvp("edge_threshold", _edge_threshold),
		cereal::make_nvp("DOG threshold ", _t)
	);
}


template<class Archive>
inline void SIFT_GPU_Image_describer::serialize(Archive & ar)
{
	ar(cereal::make_nvp("params", params_));
}

CEREAL_REGISTER_TYPE_WITH_NAME(SIFT_GPU_Image_describer, "SIFT_GPU_Image_describer");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, SIFT_GPU_Image_describer)



#endif //USE_SIFTGPU



#endif