#ifndef __POLYGON_H__
#define __POLYGON_H__

#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"

#include "DataFormat/Image2D.h"

namespace llcv {

  class Polygon {
  public:
    Polygon(){}
    Polygon(const larocv::GEO2D_Contour_t& ctor, const geo2d::Vector<float>& start);
    ~Polygon(){}
    
    int NumberDefects(float dist_thresh=0.0) const;
    int NumberDefectsNoStart(float dist_thresh=0.0) const;

    float LargestDefect() const;
    float SmallestDefect() const;

    float LargestDefectNoStart() const;
    float SmallestDefectNoStart() const;

    float EmptyAreaRatio() const;
    float EmptyArea() const;

    float PocketArea() const;
    float PocketAreaNoStart() const;

    const larocv::GEO2D_Contour_t Contour() const { return _ctor; }
    const larocv::GEO2D_Contour_t Hull() const { return _hull; }

    float Area() const;
    float Perimeter() const;
    float Charge(const larcv::Image2D& img2d, const cv::Mat& img) const;
    void DetectBranching(const cv::Mat& img,
			 const float rad,
			 const float thickness,
			 const int edge_size,
			 const int branch_size);


    const larocv::GEO2D_Contour_t Edges() const { return _edge_pt_v; }
    const larocv::GEO2D_Contour_t Branches() const { return _branch_pt_v; }

    float Fraction(const cv::Mat& img1, const cv::Mat& img2) const;

  private:
    void Construct();
    void VetoStartPoint();

  private:

    geo2d::Vector<float> _start;

    larocv::GEO2D_Contour_t _ctor;
    larocv::GEO2D_Contour_t _hull;

    std::vector<int> _hull_v;
    std::vector<cv::Vec4i> _defect_info_v;
    std::vector<float> _defect_dist_v;
    
    int _veto_ctor_id;
    int _veto_hull_id;

    larocv::GEO2D_Contour_t _edge_pt_v;
    larocv::GEO2D_Contour_t _branch_pt_v;

  };



}

#endif
