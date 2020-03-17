#ifndef UTILFUNCTIONS_H
#define UTILFUNCTIONS_H
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <math.h>

// ROOT
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"

// larcv
#include "DataFormat/EventImage2D.h"

// larutil
#include "LArUtil/LArProperties.h"
#include "LArUtil/DetectorProperties.h"
#include "LArUtil/Geometry.h"

namespace larlitecv {

  bool is_close_enough(std::vector<double> pt1, std::vector<double> compare_pt, double thresh=0.1);
  double distance_between_pt(const std::vector<double>& pt1,const std::vector<double>& pt2);
  double distance_between_pt(const TVector3& pt1 , const TVector3& pt2);
  double distance_point_plane(const TVector3& pt, double a, double b, double c, double d);
  std::vector<double> define_plane(const TVector3& pt1, const TVector3& pt2, const TVector3& pt3);
  TVector3 get_line_plane_intersection_pt(const TVector3& pt1,const TVector3& pt2, std::vector<double> planedef, bool check_between=false);
  double distance_between_pt2d(const int x1, const int y1, const int x2, const int y2);
  TVector3 get_unit_vector(const TVector3& v1, const TVector3& v2);
  std::vector<int> getProjectedPixel(const TVector3& tvec,
  				   const larcv::ImageMeta& meta,
  				   const int nplanes,
  				   const float fracpixborder=1.5 );

  std::vector<int> getProjectedPixel(const std::vector<double>& pos3d,
  				   const larcv::ImageMeta& meta,
  				   const int nplanes,
  				   const float fracpixborder=1.5 );
}
#endif
