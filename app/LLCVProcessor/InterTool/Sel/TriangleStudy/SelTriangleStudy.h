#ifndef __SELTRIANGLESTUDY_H__
#define __SELTRIANGLESTUDY_H__

#include "InterTool_Core/InterSelBase.h"
#include "LArOpenCV/ImageCluster/AlgoClass/PixelScan3D.h"
#include "TStopwatch.h"
#include "InterTool_Util/Triangle.h"
#include "InterTool_Util/Polygon.h"
#include <array>

namespace llcv {

  class SelTriangleStudy : public InterSelBase { 

  public:

  SelTriangleStudy(std::string name="SelTriangleStudy") : InterSelBase(name), _outtree(nullptr) {}
    ~SelTriangleStudy(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize();
    double Select();
    void Finalize();
    
  private:
    TTree* _outtree;

    size_t _cropx;
    size_t _cropy;
    size_t _n_neighbors;
    float _brem_dist;
    int _brem_size;
    float _min_defect_sz;

    larocv::PixelScan3D _PixelScan3D;

    TStopwatch _twatch;

  private:
    std::pair<float,float> ShowerAngle(const larlite::shower& shower);

    int DetectBrem(Triangle triangle, const larocv::GEO2D_ContourArray_t& other_ctor_v);
    
    void FillNeighbors(const std::array<cv::Mat,3>& parent_v,
		       std::array<cv::Mat,3>& child_v, 
		       const larocv::GEO2D_Contour_t& vertex_ctor) const;
    
    larocv::GEO2D_Contour_t Neighbors(const geo2d::Vector<int>& pt, const cv::Mat& mat) const;

    void ResizeOutput(size_t sz);

    void FillTriangle(const cv::Mat& img, const Triangle& tri, const size_t idx, const size_t plane);
    void FillPolygon(const Polygon& poly, const size_t idx, const size_t plane);

  private:

    std::vector<float> _triangle_height_U_v;
    std::vector<float> _triangle_height_V_v;
    std::vector<float> _triangle_height_Y_v;

    std::vector<float> _triangle_base_U_v;
    std::vector<float> _triangle_base_V_v;
    std::vector<float> _triangle_base_Y_v;

    std::vector<float> _triangle_area_U_v;
    std::vector<float> _triangle_area_V_v;
    std::vector<float> _triangle_area_Y_v;

    std::vector<float> _triangle_straight_U_v;
    std::vector<float> _triangle_straight_V_v;
    std::vector<float> _triangle_straight_Y_v;
      
    std::vector<float> _triangle_empty_area_ratio_U_v;
    std::vector<float> _triangle_empty_area_ratio_V_v;
    std::vector<float> _triangle_empty_area_ratio_Y_v;

    std::vector<float> _triangle_empty_area_U_v;
    std::vector<float> _triangle_empty_area_V_v;
    std::vector<float> _triangle_empty_area_Y_v;

    std::vector<int> _triangle_brem_U_v;
    std::vector<int> _triangle_brem_V_v;
    std::vector<int> _triangle_brem_Y_v;

    std::vector<int> _polygon_number_defects_U_v;
    std::vector<int> _polygon_number_defects_V_v;
    std::vector<int> _polygon_number_defects_Y_v;

    std::vector<int> _polygon_number_defects_no_start_U_v;
    std::vector<int> _polygon_number_defects_no_start_V_v;
    std::vector<int> _polygon_number_defects_no_start_Y_v;

    std::vector<float> _polygon_largest_defect_U_v;
    std::vector<float> _polygon_largest_defect_V_v;
    std::vector<float> _polygon_largest_defect_Y_v;
    
    std::vector<float> _polygon_smallest_defect_U_v;
    std::vector<float> _polygon_smallest_defect_V_v;
    std::vector<float> _polygon_smallest_defect_Y_v;
    
    std::vector<float> _polygon_largest_defect_no_start_U_v;
    std::vector<float> _polygon_largest_defect_no_start_V_v;
    std::vector<float> _polygon_largest_defect_no_start_Y_v;
      
    std::vector<float> _polygon_smallest_defect_no_start_U_v;
    std::vector<float> _polygon_smallest_defect_no_start_V_v;
    std::vector<float> _polygon_smallest_defect_no_start_Y_v;

    std::vector<float> _polygon_empty_area_ratio_U_v;
    std::vector<float> _polygon_empty_area_ratio_V_v;
    std::vector<float> _polygon_empty_area_ratio_Y_v;
    
    std::vector<float> _polygon_empty_area_U_v;
    std::vector<float> _polygon_empty_area_V_v;
    std::vector<float> _polygon_empty_area_Y_v;

    std::vector<float> _polygon_pocket_area_U_v;
    std::vector<float> _polygon_pocket_area_V_v;
    std::vector<float> _polygon_pocket_area_Y_v;

    std::vector<float> _polygon_pocket_area_no_start_U_v;
    std::vector<float> _polygon_pocket_area_no_start_V_v;
    std::vector<float> _polygon_pocket_area_no_start_Y_v;

  };

}


#endif
