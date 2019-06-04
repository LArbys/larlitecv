#ifndef __LINE_EXTENSION_H__
#define __LINE_EXTENSION_H__

#include "LLCVBase/llcv_base.h"

#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"

namespace llcv {

  class LineExtension : public llcv_base {

  public:
    
  LineExtension() : llcv_base("LineExtension") {}
    ~LineExtension() {} 
    
    const cv::Mat& Image() const { return _img; }
    void SetImageDimension(const cv::Mat& img, const cv::Mat& dead);
    void SetCosmicPixels(const std::vector<larocv::GEO2D_Contour_t>& cosmic_ctor_v);
    
    bool AtBoundary(const geo2d::Vector<float>& pt) const;
    bool AtParticle(const geo2d::Vector<float>& pt) const;

    bool InBoundary(const geo2d::Vector<float>& pt) const;
    bool InParticle(const geo2d::Vector<float>& pt) const;

    geo2d::Vector<float> ExtendAcross(const geo2d::Vector<float> start, const geo2d::Vector<float> edge) const;

  private:

    cv::Mat _img;
    cv::Mat _white_img;

  };

}

#endif
