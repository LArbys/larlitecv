#ifndef __TAGGER_CROI_ALGO_CONFIG_H__
#define __TAGGER_CROI_ALGO_CONFIG_H__

// larcv
#include "Base/PSet.h"

// larlitecv
#include "ThruMu/BoundaryMuonTaggerAlgoConfig.h"
#include "ThruMu/FlashMuonTaggerAlgoConfig.h"
#include "ThruMu/ThruMuTrackerConfig.h"
#include "StopMu/StopMuFilterSpacePoints.h"
#include "StopMu/StopMuClusterConfig.h"
#include "StopMu/StopMuFoxTrotConfig.h"
#include "UntaggedClustering/ClusterGroupAlgo.h"
#include "ContainedROI/TaggerFlashMatchAlgoConfig.h"

namespace larlitecv {

  class TaggerCROIAlgoConfig {
  public:
    TaggerCROIAlgoConfig();
    virtual ~TaggerCROIAlgoConfig() {};

    BoundaryMuonTaggerAlgoConfig  sidetagger_cfg;
    FlashMuonTaggerAlgoConfig     flashtagger_cfg;
    ThruMuTrackerConfig           thrumu_tracker_cfg;
    StopMuFilterSpacePointsConfig stopmu_filterpts_cfg;
    StopMuClusterConfig           stopmu_cluster_cfg;
    StopMuFoxTrotConfig           stopmu_foxtrot_cfg;
    ClusterGroupAlgoConfig        untagged_cluster_cfg;
    TaggerFlashMatchAlgoConfig    croi_selection_cfg;

    larcv::PSet input_write_cfg;
    larcv::PSet thrumu_write_cfg;
    larcv::PSet stopmu_write_cfg;
    larcv::PSet croi_write_cfg;

    bool run_thrumu_tracker;
    bool recluster_stop_and_thru;

    static TaggerCROIAlgoConfig makeConfigFromFile( std::string );

  };

}

#endif
