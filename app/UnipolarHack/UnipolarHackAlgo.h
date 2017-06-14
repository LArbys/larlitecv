#ifndef __UNIPOLARHACK_H__
#define  __UNIPOLARHACK_H__

#include <vector>
#include <deque>

#include "DataFormat/Image2D.h"

namespace larlitecv {

  class UnipolarROI_t {
  public:
    UnipolarROI_t() {
      reset();
    };
    virtual ~UnipolarROI_t() {};
    
    int neg_peak;
    int pos_peak;    
    int start_row;
    int end_row;
    int plat_row;
    float plat_val;
    int nplat;
    std::vector<float> plat_vals;
    float neg_peak_val;
    float pos_peak_val;
    void reset() {
      neg_peak = -1;
      pos_peak = -1;      
      start_row = -1;
      end_row = -1;
      plat_row = -1;
      plat_val = 0;
      nplat = 0;
      plat_vals.clear();
      plat_vals.reserve(200);
      neg_peak_val = 0.;
      pos_peak_val = 0.;
    };
  };
  
  class UnipolarHackAlgo {
  public:
    UnipolarHackAlgo() {};
    virtual ~UnipolarHackAlgo() {};


    std::vector< larcv::Image2D > hackUnipolarArtifact( const std::vector<larcv::Image2D>& img_v, const std::vector<int>& applytoplane, const std::vector<float>& neg_threshold );

    void scanForPulse( const int col, const larcv::Image2D& img, const float neg_threshold, std::vector<UnipolarROI_t>& rois );

  };
  
}

#endif
