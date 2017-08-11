#ifndef __BMTCombinedFilter_h__
#define __BMTCombinedFilter_h__

/* ------------------------------------------------------
 * BMTCombinedFilter
 * 
 * This class runs the different end point filters.
 * Allows me to develop the filters while keeping the main
 * tagger code unaffected.
 *
 * author: Taritree (twongj01@tufts.edu)
 *
 * history:
 * 2017/08/03 - first writing
 * 
 * -----------------------------------------------------*/

#include "BMTCombinedFilterConfig.h"
#include "TaggerContourTools/BMTCV.h"

namespace larlitecv {

  class BMTCombinedFilter {

  public:
    
    BMTCombinedFilter() {}; // dont use this. but needed for dictionary
    BMTCombinedFilter( const larlitecv::BMTCombinedFilterConfig& config );
    virtual ~BMTCombinedFilter() {};

    void filterBMTEndpoints( const std::vector< larlitecv::BoundarySpacePoint* >& endpt_v,
			     const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
			     const larlitecv::BMTCV& contours, std::vector<int>& decision_v );

    const larlitecv::BMTCombinedFilterConfig& getConfig() const { return m_config; };
    
  protected:
    
    larlitecv::BMTCombinedFilterConfig m_config;


  };

}

#endif
