#ifndef __LINEFOLLOW_H__
#define __LINEFOLLOW_H__

#include "LLCVBase/llcv_base.h"

#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"

#include "LArOpenCV/ImageCluster/AlgoClass/DeadWirePatch.h"
#include "LArOpenCV/ImageCluster/AlgoClass/PiRange.h"

namespace llcv {

  class LineFollow : public llcv_base {

  public:
    LineFollow();
    ~LineFollow() {}

    larocv::GEO2D_ContourArray_t FollowEdgeLine(const geo2d::Vector<float>& start);

    void SetImageDimension(const cv::Mat& img, const cv::Mat& dead);

    larocv::GEO2D_Contour_t EdgePoints();

    cv::Mat& Image()      { return _img;       }
    cv::Mat& DeadImage()  { return _dead_img;  }
    cv::Mat& BondImage()  { return _bond_img;  }
    cv::Mat& WhiteImage() { return _white_img; }
    cv::Mat& BlackImage() { return _black_img; }

  private:

    cv::Mat _img;
    cv::Mat _dead_img;
    cv::Mat _bond_img;
    cv::Mat _black_img;
    cv::Mat _white_img;
    
    int _thickness;
    float _radius;
    std::vector<float> _radius_v;

    larocv::PiRange _PiRange;
    larocv::DeadWirePatch _DeadWirePatch;

  private:

    bool InitializeFirstPoint(geo2d::Vector<float> start, geo2d::Vector<float>& init_pt);
    larocv::GEO2D_Contour_t AsLine(geo2d::Vector<float> pt1, geo2d::Vector<float> pt2);
    float DistanceToEdge(const geo2d::Vector<float>& pt) const;

    int Quadrant(const geo2d::Vector<float>& pt) const;
    float Angle(const geo2d::Vector<float>& origin, const geo2d::Vector<float>& pt1) const;

  private:

    
  };

}

#endif
