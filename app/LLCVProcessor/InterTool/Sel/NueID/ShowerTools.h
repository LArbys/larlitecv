#ifndef __SHOWERTOOLS_H__
#define __SHOWERTOOLS_H__

#include "LLCVBase/llcv_base.h"

#include "Object2D.h"
#include "DataFormat/Image2D.h"
#include "InterTool_Util/TruncMean.h"

#include <array>

namespace llcv {
  
  class ShowerTools : public llcv_base {
  
  public:
    
  ShowerTools() : llcv_base("ShowerTools") {}
    ~ShowerTools() {}
    
    void ReconstructAngle(const std::vector<larcv::Image2D*>& img_v,
			  const std::array<cv::Mat,3>& aimg_v, 
			  Object2DCollection& obj_col);
    
    void ReconstructLength(const std::vector<larcv::Image2D*>& img_v,
			   const std::array<cv::Mat,3>& aimg_v,
			   Object2DCollection& obj_col);

    void ReconstructdQdx(const std::vector<larcv::Image2D*>& img_v,
			 const std::array<cv::Mat,3>& aimg_v,
			 Object2DCollection& obj_col,
			 double dtrunk);
    
    void ReconstructdQdxProfile(const std::vector<larcv::Image2D*>& img_v,
				const std::array<cv::Mat,3>& aimg_v,
				Object2DCollection& obj_col);

    void TruncatedQdxProfile(Object2DCollection& obj_col, const float ftsigma);
    
    std::array<float,3> ComputePCA(std::vector<std::array<float,3> > pts_v, 
				   const Object2DCollection& obj_col);


    TruncMean _TruncMean;    

  };


}

#endif
