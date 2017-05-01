#ifndef __STOPMU_FOXTROT_H__
#define __STOPMU_FOXTROT_H__

/* =============================================================
 *
 * StopMuFoxTrot
 *
 * Bridge between tagger code and the FoxTrotTrackerAlgo
 *
 * ============================================================*/

#include <vector>

// larcv
#include "DataFormat/Image2D.h"

#include "TaggerTypes/BMTrackCluster3D.h"
#include "TaggerTypes/BoundarySpacePoint.h"
#include "ChargeSegmentAlgos/FoxTrotTrackerAlgo.h"

#include "StopMuFoxTrotConfig.h"


namespace larlitecv {

  class StopMuFoxTrot {
  public:
    StopMuFoxTrot( const StopMuFoxTrotConfig& config );
    virtual ~StopMuFoxTrot();

    std::vector< larlitecv::BMTrackCluster3D > findStopMuTracks( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
								 const std::vector<larcv::Image2D>& thrumu_v, const std::vector<BoundarySpacePoint>& endpts_v );

  protected:
    
    StopMuFoxTrotConfig m_config;
    FoxTrotTrackerAlgo* m_algo;
    
  };

}

#endif
