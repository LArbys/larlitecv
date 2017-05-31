#ifndef __TAGGER_CROI_ALGO_H__
#define __TAGGER_CROI_ALGO_H__

#include "TaggerCROIAlgoConfig.h"
#include "TaggerCROITypes.h"

namespace larlitecv {

  class TaggerCROIAlgo {

    TaggerCROIAlgo() {};

  public:

    TaggerCROIAlgo( const TaggerCROIAlgoConfig& config );
    virtual ~TaggerCROIAlgo() {};


    ThruMuPayload runThruMu( const InputPayload& data );
    StopMuPayload runStopMu( const InputPayload& input, const ThruMuPayload& thrumu );
    CROIPayload   runCROISelection( const InputPayload& input, const ThruMuPayload& thrumu, const StopMuPayload& stopmu );    

    void printTimeTracker( int num_events=0 );

    // subroutines: for some development and flexibility
    void runBoundaryTagger( const InputPayload& data, ThruMuPayload& thrumu );
    void runThruMuTracker( const InputPayload& data, ThruMuPayload& thrumu );
    
  protected:

    TaggerCROIAlgoConfig m_config;
    std::vector<float> m_time_tracker;
    enum Stages_t { kThruMuConfig = 0, kThruMuBMT, kThruMuFlash, kThruMuTracker, kStopMuTracker, kUntagged, kCROI, kNumStages };

  };
  

}

#endif
