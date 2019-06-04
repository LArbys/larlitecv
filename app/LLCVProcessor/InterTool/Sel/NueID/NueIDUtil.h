#ifndef __NUEIDUTIL_H__
#define __NUEIDUTIL_H__

#include "TStopwatch.h"
#include "ContourScan.h"
#include "InterTool_Util/Triangle.h"
#include "MatchObjectAlgoTimeIOU.h"
#include "Object2D.h"
#include "CosmicTag.h"
#include "LineExtension.h"
#include "ShowerTools.h"
#include "LArOpenCV/ImageCluster/AlgoClass/PixelScan3D.h"
#include "DataFormat/track.h"

namespace llcv {

  size_t FindClosestContour(const larocv::GEO2D_ContourArray_t& ctor_v,
			    const geo2d::Vector<float>& pt,
			    float& distance);
  
  larocv::GEO2D_Contour_t MaximizeLine(const cv::Mat& timg3d_par,
				       const Triangle& triangle,
				       float& nline_pixels,
				       geo2d::Vector<float>& edge,
				       cv::Mat& white_img);
  
  larocv::GEO2D_Contour_t MaximizeTriangleLine(const cv::Mat& timg3d_mask,
					       const Triangle& triangle,
					       float& nline_pixels,
					       float& npar_pixels,
					       geo2d::Vector<float>& edge,
					       cv::Mat& white_img);
  
  larocv::GEO2D_Contour_t MaximizePolygonLine(const cv::Mat& timg3d_mask,
					      const std::vector<Polygon>& polygon_v,
					      Triangle& triangle,
					      float& nline_pixels,
					      float& npar_pixels,
					      geo2d::Vector<float>& edge,
					      cv::Mat& white_img);

  float MinimizeToEdge(const larocv::GEO2D_Contour_t& ctor,
		       size_t crop_x,
		       size_t crop_y);
    
  float NearestPolygonToCosmic(const std::vector<Polygon>& polygon_v,
			       const std::vector<CosmicTag>& CosmicTag_v,
			       size_t plane);
  
  float NearestPolygonToCosmicEnd(const std::vector<Polygon>& polygon_v,
				  const std::vector<CosmicTag>& CosmicTag_v,
				  size_t plane);
  
  float PointCosmicDistance(const geo2d::Vector<float>& pt,
			    const std::vector<CosmicTag>& CosmicTag_v,
			    const size_t plane);
  
  float PointCosmicEndDistance(const geo2d::Vector<float>& pt,
			       const std::vector<CosmicTag>& CosmicTag_v,
			       const size_t plane);

  int DetectBrem(Triangle& triangle, 
		 const larocv::GEO2D_ContourArray_t& other_ctor_v,
		 std::vector<size_t>& id_v,
		 std::vector<size_t>& inside_v,
		 float brem_dist,
		 float brem_size);
    
  std::array<float,3> EstimateDirection(const Object2DCollection& obj_col,
					cv::Mat& white_img,
					ContourScan& ContourScan_,
					ShowerTools& ShowerTools_);

  void ComputeLineParameters(Object2DCollection& obj_col,
			     const std::array<cv::Mat,3>& cimg_v,
			     cv::Mat& white_img);
  
  void SplitLineParameters(Object2DCollection& obj_col,
			   const std::array<cv::Mat,3>& cimg_v,
			   cv::Mat& white_img);

  void FillTrack(const Object2DCollection& obj_col,
		 larlite::track& out_track);

}

#endif
