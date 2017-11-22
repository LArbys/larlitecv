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

// larlite
#include "FhiclLite/FhiclLiteUtilFunc.h"

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

    int verbosity;
    bool run_thrumu_tracker;
    bool use_truth_endpoints;
    bool recluster_stop_and_thru;

    std::string larcv_image_producer;
    std::string larcv_chstatus_producer;
    bool DeJebWires;
    float jebwiresfactor;
    std::vector<float> emptych_thresh;
    std::string chstatus_datatype;
    std::vector<std::string> opflash_producers;
    bool RunThruMu;
    bool RunStopMu;
    bool RunCROI;
    bool save_thrumu_space;
    bool save_stopmu_space;
    bool save_croi_space;
    bool save_mc;
    bool load_mctrack;
    std::string mctrack_producer;
    std::string trigger_producer;    
    bool skip_empty_events;
    bool apply_unipolar_hack;

    static TaggerCROIAlgoConfig makeConfigFromFile( std::string );

  };

}

#endif
