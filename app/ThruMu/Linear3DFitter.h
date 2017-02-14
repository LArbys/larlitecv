#ifndef __LINEAR_3D_FITTER_H__
#define __LINEAR_3D_FITTER_H__

#include <vector>

// LArCV
#include "DataFormat/Image2D.h"
#include "Base/PSet.h"

// LArLiteCV
#include "BMTrackCluster3D.h"
#include "BoundarySpacePoint.h"

namespace larlitecv {

  class PointInfo {
  public:
    std::vector<float> xyz;
    std::vector<int> wire_id;
    int row;
    float tick;
    std::vector<int> cols;
    std::vector<bool> planehascharge;
    std::vector<bool> planehasbadch;
    bool goodpoint;
    int planeswithcharge;
    PointInfo() {
      planeswithcharge = 0;			
      goodpoint = false;
      tick = 0.0;
      row = 0;
    };
    virtual ~PointInfo() {};

    bool operator==(const PointInfo& rhs) const {
      if ( row!=rhs.row ) return false;
      if ( cols.size()!=rhs.cols.size() ) return false;
      for ( int i=0; i<(int)cols.size(); i++ ) {
        if ( cols[i]!=rhs.cols[i] ) return false;
      };
      return true;
    };
  };

  class PointInfoList : public std::vector<PointInfo> {
  public:
    PointInfoList() {
      num_pts_w_allcharge = 0;
      num_pts_w_majcharge = 0;
      num_pts_w_allbadch = 0;
      num_pts_w_allempty = 0;
      num_pts_good = 0;
    };
    virtual ~PointInfoList() {};

#ifndef __CINT__
#ifndef __CLING__
    // we must hide this from ROOT's interpretors
    void emplace( PointInfo&& pt );
#endif
#endif

    float fractionHasChargeWith3Planes() { return float(num_pts_w_allcharge)/float(size()); };
    float fractionHasBadChOn3Planes() { return float(num_pts_w_allbadch)/float(size()); };
    float fractionHasNoChargeOn3Planes() { return float(num_pts_w_allempty)/float(size()); };
    float fractionHasChargeOnMajorityOfPlanes() { return float(num_pts_w_majcharge)/float(size()); };
    float fractionGood() { return float(num_pts_good)/float(size()); };

    int num_pts_w_allcharge;
    int num_pts_w_majcharge;
    int num_pts_w_allbadch;
    int num_pts_w_allempty;
    int num_pts_good;
  };

  class Linear3DFitterConfig {
  public:
    Linear3DFitterConfig();
    virtual ~Linear3DFitterConfig() {};

    float trigger_tpc_tick;
    float min_ADC_value;
    float step_size;
    int neighborhood_square;
    int neighborhood_posttick;// allows us to extend in later ticks to account for space-charge effect delay

    static Linear3DFitterConfig makeFromPSet( const larcv::PSet& pset );
    
  };


  class Linear3DFitter  {
  public:
    Linear3DFitter( const Linear3DFitterConfig& cfg ) : m_config(cfg) {}; 
    virtual ~Linear3DFitter() {};

    Linear3DFitterConfig m_config;

    PointInfoList findpath( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const int start_row, const int goal_row, const std::vector<int>& start_cols, const std::vector<int>& goal_cols  );

    PointInfoList pointsOnTrack(const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const std::vector<float>& initial_point_coords, const std::vector<float>& final_point_coords, 
     const float step_size, const float min_ADC_value, const int neighborhood_size);

    bool doesNeighboringPixelHaveCharge(const larcv::Image2D& img, int central_row, int central_col, int neighborhood_size, float min_ADC_value);

    bool isNeighboringPixelDeadPixel(const larcv::Image2D& badch_v, int central_row, int central_col, int neighborhood_size);

    BMTrackCluster3D makeTrackCluster3D( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const BoundarySpacePoint& start_pt, const BoundarySpacePoint& end_pt, const PointInfoList& infolist );

    void getTrackExtension( const PointInfoList& infolist, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
      const float max_extension_length, PointInfoList& start_extension, PointInfoList& end_extension );


  };


}

#endif