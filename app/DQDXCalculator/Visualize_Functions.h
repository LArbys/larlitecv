#ifndef VISUALIZE_FUNCTIONS_H
#define VISUALIZE_FUNCTIONS_H
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

// larlite
#include "DataFormat/track.h"
// larutil
#include "LArUtil/Geometry.h"
//DQDXCalculator
#include "UtilFunctions.h"

namespace larlitecv {
  void make_evdisp(TH2D& hist, std::string title, std::string x_title="Col (Wire)", std::string y_title="Row (6 Ticks)", double maxval=-1, TH2D* vtx_h=NULL);
  void make_evdisp_single(const larcv::Image2D& img, std::string title, std::string x_title="Col (Wire)", std::string y_title="Row (6 Ticks)", double maxval =-1, TH2D* vtx_h=NULL);
  void make_evdisp_triple(const std::vector<larcv::Image2D>& img_v,  std::string title, std::string x_title="Col (Wire)", std::string y_title="Row (6 Ticks)", double maxval=-1, TH2D* vtx_h=NULL);
  void make_dqdx_curve(const larlite::track& track, int plane, std::string filename="dqdx_curve");
  void make_dqdx_curve(const std::vector<double>& dqdx_v, const std::vector<double>& distance_v, std::string filename="dqdx_curve");
  void print_signal();
  void print_tvec(const TVector3 tvec);
}
#endif
