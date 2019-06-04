#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"

namespace llcv {

  class Triangle {
  public:
    Triangle(){}
    Triangle(const larocv::GEO2D_Contour_t& ctor,
	     const geo2d::Vector<float>& start);
    Triangle(const geo2d::Vector<float>& apex, 
	     const geo2d::Vector<float>& base1,
	     const geo2d::Vector<float>& base2) 
      : _apex(apex), _base_pt1(base1), _base_pt2(base2) {}

    ~Triangle(){}
    
    void Expand(const cv::Mat& img, const float fraction);
    void Tighten(const cv::Mat& img, const float radius, const float fraction);
    float StraightLineTest(const cv::Mat& img) const;
    void Extend(const float fraction);
    Triangle RotateToPoint(const geo2d::Vector<float>& pt,float scale) const;
    void SetContour(const larocv::GEO2D_Contour_t& ctor) { _ctor = ctor; }


  private:
    void Construct();
    void MovePt(const cv::Mat& ctor_img, geo2d::Vector<float>& base_pt);
    geo2d::Vector<float> MidPoint(const geo2d::Vector<float>& pt1,const geo2d::Vector<float>& pt2) const;

  public:
    const larocv::GEO2D_Contour_t& Contour() const { return _ctor; }
    larocv::GEO2D_Contour_t AsContour() const;
    const geo2d::Vector<float>& Apex() const { return _apex;}
    const geo2d::Vector<float>& Base1() const { return _base_pt1; }
    const geo2d::Vector<float>& Base2() const { return _base_pt2; }
    geo2d::Vector<float> MidPoint() const { return MidPoint(_base_pt1,_base_pt2); }
    
    float Height() const;

    float EmptyAreaRatio() const;
    float EmptyArea() const;
    
    float BaseLength() const;
    float Area() const;

    float Coverage(const cv::Mat& img) const;
    

  private:
    geo2d::Vector<float> _apex;
    geo2d::Vector<float> _base_pt1;
    geo2d::Vector<float> _base_pt2;

    larocv::GEO2D_Contour_t _ctor;

  };

}

#endif
