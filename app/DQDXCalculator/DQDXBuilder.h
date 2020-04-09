#ifndef DQDXBUILDER_H
#define DQDXBUILDER_H

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <math.h>

// larutil
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

// larcv
#include "DataFormat/EventImage2D.h"

// larlite
#include "DataFormat/track.h"

//DQDXCalculator Package
#include "DQDXCalculator/UtilFunctions.h"
#include "DQDXCalculator/Visualize_Functions.h"
#include "LArUtil/DetectorProperties.h"

namespace larlitecv {

  /**
      \class DQDXBuilder
      Create this class to feed larlite reco3d tracks and have it return
      larlite tracks with dqdx calculations
   */


   class DQDXBuilder {
   public:
     DQDXBuilder();
     DQDXBuilder(std::vector<larcv::Image2D> img_v);

     ~DQDXBuilder() {}



     larlite::geo::View_t views[3] ;
     std::vector<std::vector<std::vector<double> > >  wire_planar_equations; // Plane < Wire < a,b,c,d > > Where abcd are plane equation constants
     const larutil::Geometry* geo ;
     TH2D residual_dqdx;
     std::vector<TVector3> new_steps_y;
     std::vector<double> dq_y;
     std::vector<double> dx_y;
     std::vector<double> s_distance_y;
     std::vector<larcv::Image2D> get_img_v() {return _img_v;} //Returns copy of _img_v
     void Initialize_DQDXBuilder();
     void set_img_v(std::vector<larcv::Image2D> img_v) {_img_v = img_v;} // Sets _img_v
     double sum_charge_dq(const int low_row, const int high_row, const int col, const int plane=2, const int buffer=4);
     void Initialize_Wirept_Vector();
     int calc_row_from_x(double x1, const larcv::ImageMeta& meta);
     // larlite::track calc_dqdx_track(const larlite::track& reco3d_track);
     void Add_Track_To_ResidualDQDX(const larlite::track& reco3d_track, const int plane);
     void Plot_Orig_New_Raw_Track_Crop(const larlite::track& orig_track, const larlite::track& new_track);

     larlite::track calc_dqdx_track_revamp(const larlite::track& reco3d_track, const int plane);
     TVector3 get_halfway_plane_pt( int w1,  int w2, const TVector3 pt1, const TVector3 pt2, const int plane);
     void test_wirecoord_driftplanes(larlite::track reco3d_track);

   private:
     std::vector<larcv::Image2D> _img_v; // The three plane views in adc


   };
   //End of Class


}
#endif
