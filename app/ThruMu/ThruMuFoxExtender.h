#ifndef __THRUMUFOXEXTENDER_H__
#define __THRUMUFOXEXTENDER_H__

/* ======================================================
 *  FoxTrotTrackerAlgoConfig
 *
 *  Wrapper for FoxTrotTracker with an aim towards
 *  extending thrumu tracks. Might not be necessary.
 *
 * ======================================================*/

#include <vector>

#include "DataFormat/Image2D.h"

#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "ChargeSegmentAlgos/FoxTrotTrackerAlgoConfig.h"
#include "ChargeSegmentAlgos/FoxTrotTrackerAlgo.h"

namespace larlitecv {

  class ThruMuFoxExtender {
  public:
    ThruMuFoxExtender( const FoxTrotTrackerAlgoConfig& config )
    : m_config(config)
    {};
    virtual ~ThruMuFoxExtender() {};

    void setConfig( const FoxTrotTrackerAlgoConfig& config ) {
    	m_config = config;
    };

    bool extendTrack( std::vector<std::vector<double> >& track, const std::vector<larcv::Image2D>& img_v,
    	const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v );


  protected:

  	enum ExtEnd_t { kBack=0, kFront };

  	FoxTrotTrackerAlgoConfig m_config;
  	FoxTrack m_extensions[2]; // back and front

    FoxTrack extendFromPoint( const std::vector<float>& pos, const std::vector<float>& dir, const std::vector<larcv::Image2D>& img_v,
    													const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v );

  	bool appendExtension( ThruMuFoxExtender::ExtEnd_t end, FoxTrack& extension,
			std::vector<std::vector<double> >& track, const float max_step_size, const std::vector<larcv::Image2D>& img_v,
    	const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
    	const std::vector<float>& thresholds, const int hit_neighborhood );


  };

}

#endif
