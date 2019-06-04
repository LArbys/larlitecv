#ifndef __OBJECT2D_H__
#define __OBJECT2D_H__

#include "InterTool_Util/Triangle.h"
#include "InterTool_Util/Polygon.h"

#include "TVector3.h"

#include "DataFormat/Image2D.h"

namespace llcv {

  class Object2D {
  public:
    Object2D(){}
    ~Object2D() {}
    

  public:
    float LineLength() const { return geo2d::dist(_triangle.Apex(),_edge); }
    const geo2d::Vector<float>& Start() const { return _triangle.Apex(); }
    const geo2d::Vector<float>& Edge() const { return _edge; }
    int NPolygons() const { return (int)(_polygon_v.size()); }
    float LineFrac() const { return _line_frac; }
    float LineFracEmpty(const cv::Mat& white_img) const;

    const Triangle& triangle() const { return _triangle; }
    const Triangle& brem_triangle() const { return _brem_triangle; }
    const std::vector<Polygon>& Polygons() const { return _polygon_v; }
    const std::vector<Polygon>& ExpandedPolygons() const { return _expand_polygon_v; }
    const larocv::GEO2D_Contour_t& Line() const { return _line; }
    const geo2d::Vector<float> Dir() const { return _dir; }

    size_t Plane() const { return _plane; }

    int NBrem() const { return _n_brem; }

    float LinedX() const;
    float LinedY() const;
    
    float Charge(const larcv::Image2D& img2d, const cv::Mat& img) const;
    
    float Length() const { return _length; }

    double dQdx() const { return _dqdx; }
    
    const std::vector<float>& dQdxProfile() const { return _dqdx_v; }
    const std::vector<float>& TdQdxProfile() const { return _tdqdx_v; }
    const std::vector<float>& dxProfile() const { return _dx_v; }

    float dQdxStep() const { return _dqdx_step; }

    float dQdxPitch() const { return _dqdx_pitch; }
    
    int BremIndex() const { return _brem_index; }
    
    float Fraction(const cv::Mat& img1, const cv::Mat& img2) const;
    
    void LineVertex(const larcv::Image2D& img2d,
		    const cv::Mat& img,
		    const cv::Mat& white_img, 
		    float radius);

    float LineVertexDensity() const { return _vtx_density; }
    float LineVertexCoverage() const { return _vtx_coverage; }
    float LineVertexCharge() const { return _vtx_charge; }

    float LineMeanDist() const { return _line_mean_dist; }
    float LineMaxDist() const { return _line_max_dist; }
    float LineFirstHalfLineFrac() const { return _line_first_half_linefrac; }
    float LineSecondHalfLineFrac() const { return _line_second_half_linefrac; }

  public:

    Triangle _triangle;
    Triangle _brem_triangle;
    
    std::vector<Polygon> _polygon_v; 
    std::vector<Polygon> _expand_polygon_v;

    larocv::GEO2D_Contour_t _line;
    geo2d::Vector<float> _edge;
    float _line_frac;
    size_t _plane;
    int _n_brem;

    float _length;
    double _dqdx;

    std::vector<float> _dqdx_v;
    std::vector<float> _tdqdx_v;
    std::vector<float> _dx_v;

    float _dqdx_step;
    float _dqdx_pitch;

    int _brem_index;
    
    geo2d::Vector<float> _dir;

    float _vtx_density;
    float _vtx_coverage;
    float _vtx_charge;

    larocv::GEO2D_Contour_t _vtx_pt_v;

    float _line_mean_dist;
    float _line_max_dist;
    
    float _line_first_half_linefrac;
    float _line_second_half_linefrac;

  };

  class Object2DCollection : public std::vector<Object2D> {

  public:
    Object2DCollection() {}
    ~Object2DCollection() {}

    bool HasObject(size_t plane) const;
    const Object2D& PlaneObject(size_t plane) const;
    Object2D& PlaneObjectRW(size_t plane);
    
    void SetTheta(float theta) { _theta = theta; }
    void SetPhi(float phi)     { _phi = phi;     }

    void SetStart(float x, float y, float z)  { _start = TVector3(x,y,z); }

    void SetLength(float length) { _length = length; };
    void SetScore(float score) { _score = score; }

    const TVector3& Start() const { return _start; }
    float Theta() const { return _theta; }
    float Phi()   const { return _phi; }
    float Length() const { return _length; }
    float Score() const { return _score; }

    float dX() const { return std::cos(_theta) * std::sin(_phi); }
    float dY() const { return std::sin(_theta);                  }
    float dZ() const { return std::cos(_theta) * std::cos(_phi); }
    
    std::vector<int> Planes() const;

    std::vector<int> XDead(const std::array<cv::Mat,3>& dimg_v,
			   const cv::Mat& white_img,
			   float radius=3) const;
    
    void SetddX(float x) { _ddx = x; }
    void SetddY(float y) { _ddy = y; }
    void SetddZ(float z) { _ddz = z; }
    
    float ddX() const { return _ddx; }
    float ddY() const { return _ddy; }
    float ddZ() const { return _ddz; }
    
    std::vector<float> EndPoint() const;
    
    void SetEndX(float ex) { _end_x = ex; }
    void SetEndY(float ey) { _end_y = ey; }
    void SetEndZ(float ez) { _end_z = ez; }

  private:

    float _theta;
    float _phi;
    float _length;
    TVector3 _start;
    float _score;

    float _ddx;
    float _ddy;
    float _ddz;

    float _end_x;
    float _end_y;
    float _end_z;

  };

}

#endif
