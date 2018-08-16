#ifndef __TAGGER_CROI_ALGO_H__
#define __TAGGER_CROI_ALGO_H__

// larlitecv
#include "TaggerCROIAlgoConfig.h"
#include "TaggerCROITypes.h"
#include "TaggerContourTools/BMTCV.h"
#include "Base/DataCoordinator.h"
#include "MCTruthTools/crossingPointsAnaMethods.h"
#include "SCE/SpaceChargeMicroBooNE.h"

namespace larlitecv {

  class TaggerCROIAlgo {

    TaggerCROIAlgo() {};

  public:

    TaggerCROIAlgo( const TaggerCROIAlgoConfig& config );
    virtual ~TaggerCROIAlgo() {};


    InputPayload  loadInput( DataCoordinator& dataco );
    ThruMuPayload runBoundaryPointFinder( const InputPayload& data );
    void runThruMu( const InputPayload& data, ThruMuPayload& output );
    StopMuPayload runStopMu( const InputPayload& input, const ThruMuPayload& thrumu );
    CROIPayload   runCROISelection( const InputPayload& input, const ThruMuPayload& thrumu, const StopMuPayload& stopmu );    

    void printTimeTracker( int num_events=0 );

    // subroutines: for some development and flexibility
    void runBoundaryTagger( const InputPayload& data, ThruMuPayload& thrumu );
    void runThruMuTracker( const InputPayload& data, ThruMuPayload& thrumu );

    // mc-truth subroutines
    void runTruthBoundaryTagger( const InputPayload& input, ThruMuPayload& output );
    void runTruthThruMu( const InputPayload& input, ThruMuPayload& output );

  protected:


    TaggerCROIAlgoConfig m_config;
    std::vector<float> m_time_tracker;
    enum Stages_t { kThruMuConfig = 0, kThruMuContour, kThruMuBMT, kThruMuFlash, kThruMuFilter, kThruMuTracker, kStopMuTracker, kRecluster, kUntagged, kPCAmerge, kCROI, kNumStages };
    
    // Contour tool
    larlitecv::BMTCV m_bmtcv_algo;

    // Truth data (if used)
    CrossingPointAnaData_t m_truthxingdata;

    SpaceChargeMicroBooNE m_sce;

  };
  

}

#endif
