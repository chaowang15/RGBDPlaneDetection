#ifndef PLANEDETECTION_H
#define PLANEDETECTION_H

#include <iostream>
#include "opencv2/opencv.hpp"
#include <string>
#include <fstream>
#include <Eigen/Eigen>
#include "AHCPlaneFitter.hpp"
#include "mrf.h"
#include "GCoptimization.h"

using namespace std;

typedef Eigen::Vector3d VertexType;
typedef Eigen::Vector2i PixelPos;

const int kNeighborRange = 5; // boundary pixels' neighbor range
const int kScaleFactor = 5; // scale coordinate unit in mm
const float kInfVal = 1000000; // an infinite large value used in MRF optimization

// Camera intrinsic parameters.
// All BundleFusion data uses the following parameters.
const double kFx = 583;
const double kFy = 583;
const double kCx = 320;
const double kCy = 240;
const int kDepthWidth = 640;
const int kDepthHeight = 480;

struct ImagePointCloud
{
	vector<VertexType> vertices; // 3D vertices
	int w, h;

	inline int width() const { return w; }
	inline int height() const { return h; }
	inline bool get(const int row, const int col, double &x, double &y, double &z) const {
		const int pixIdx = row * w + col;
		z = vertices[pixIdx][2];
		// Remove points with 0 or invalid depth in case they are detected as a plane
		if (z == 0 || _isnan(z)) return false;
		x = vertices[pixIdx][0];
		y = vertices[pixIdx][1];
		return true;
	}
};

// Data for sum of points on a same plane
struct SumStats
{
	double sx, sy, sz; // sum of x/y/z
	double sxx, syy, szz, sxy, syz, sxz; // sum of x*x/y*y/z*z/x*y/y*z/x*z
	SumStats(){
		sx = sy = sz = sxx = syy = szz = sxy = syz = sxz = 0;
	}
};

class PlaneDetection
{
public:
	ImagePointCloud cloud;
	ahc::PlaneFitter< ImagePointCloud > plane_filter;
	vector<vector<int>> plane_vertices_; // vertex indices each plane contains
	cv::Mat seg_img_; // segmentation image
	cv::Mat color_img_; // input color image
	int plane_num_;

	/* For MRF optimization */
	cv::Mat opt_seg_img_;
	cv::Mat opt_membership_img_; // optimized membership image (plane index each pixel belongs to)
	vector<bool> pixel_boundary_flags_; // pixel is a plane boundary pixel or not
	vector<int> pixel_grayval_;
	vector<cv::Vec3b> plane_colors_;
	vector<SumStats> sum_stats_, opt_sum_stats_; // parameters of sum of points from the same plane
	vector<int> plane_pixel_nums_, opt_plane_pixel_nums_; // number of pixels each plane has

public:
	PlaneDetection();
	~PlaneDetection();

	//bool readIntrinsicParameterFile(string filename);

	bool readColorImage(string filename);

	bool readDepthImage(string filename);

	bool runPlaneDetection();

	void prepareForMRF();

	void writeOutputFiles(string filename_prefix, bool run_mrf = false);

	void writePlaneData(string filename, bool run_mrf = false);

	/************************************************************************/
	/* For MRF optimization */
	inline MRF::CostVal dCost(int pix, int label)
	{
		return pixel_boundary_flags_[pix] ? 1 : (label == plane_filter.membershipImg.at<int>(pix / kDepthWidth, pix % kDepthWidth) ? 1 : kInfVal);
	}
	inline MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
	{
		int gray1 = pixel_grayval_[pix1], gray2 = pixel_grayval_[pix2];
		return i == j ? 0 : exp(-MRF::CostVal(gray1 - gray2) * (gray1 - gray2) / 900); // 900 = sigma^2 by default
	}
	/************************************************************************/

private:
	inline int RGB2Gray(int x, int y)
	{
		return int(0.299 * color_img_.at<cv::Vec3b>(x, y)[2] +
			0.587 * color_img_.at<cv::Vec3b>(x, y)[1] +
			0.114 * color_img_.at<cv::Vec3b>(x, y)[0] +
			0.5);
	}

	void computePlaneSumStats(bool run_mrf = false);

};


#endif