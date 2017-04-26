#ifndef __FOX_TROT_TRACKER_ALGO_H__
#define __FOX_TROT_TRACKER_ALGO_H__

/* ==================================================
 *
 * FoxTrotTrackerAlgo
 *
 * Tracker using box steps.
 *
 ==================================================*/

#include <vector>

//larcv
#include "DataFormat/Image2D.h"

// larlitecv
#include "ThruMu/BoundarySpacePoint.h"
#include "ThruMu/BMTrackCluster3D.h"
#include "ChargeSegmentAlgos/RadialSegmentSearch.h"

#include "FoxTrotTrackerAlgoConfig.h"
#include "FoxTrotTrackerAlgoTypes.h"

namespace larlitecv {

  class FoxTrotTrackerAlgo {

    FoxTrotTrackerAlgo() {};

  public:

    FoxTrotTrackerAlgo( const FoxTrotTrackerAlgoConfig& config );
    virtual ~FoxTrotTrackerAlgo() {};

    FoxTrack followTrack( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
			  const BoundarySpacePoint& start );

  protected:

    FoxTrotTrackerAlgoConfig m_config;

    RadialSegmentSearch m_radialalgo;

    FoxStep getNextStep( const FoxStep& current, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, float min_dcos );

  };

}

#endif
