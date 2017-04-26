#include "TaggerCROIAlgoConfig.h"

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

namespace larlitecv {

  TaggerCROIAlgoConfig::TaggerCROIAlgoConfig()
    : input_write_cfg(larcv::PSet("input")),
      thrumu_write_cfg(larcv::PSet("thrumu")),
      stopmu_write_cfg(larcv::PSet("stopmu")),
      croi_write_cfg(larcv::PSet("croi"))
  {
    run_thrumu_tracker = true;
  }

  TaggerCROIAlgoConfig TaggerCROIAlgoConfig::makeConfigFromFile( std::string filepath ) {

    TaggerCROIAlgoConfig cfg;

    larcv::PSet root_pset           = larcv::CreatePSetFromFile( filepath );
    larcv::PSet tagger_pset         = root_pset.get<larcv::PSet>("TaggerCROI");
    larcv::PSet sidetagger_pset     = tagger_pset.get<larcv::PSet>("BMTSideTagger");
    larcv::PSet flashtagger_pset    = tagger_pset.get<larcv::PSet>("BMTFlashTagger");
    larcv::PSet thrumu_tracker_pset = tagger_pset.get<larcv::PSet>("ThruMuTracker");
    larcv::PSet stopmu_filter_pset  = tagger_pset.get<larcv::PSet>("StopMuSpacePointsFilter");
    larcv::PSet stopmu_cluster_pset = tagger_pset.get<larcv::PSet>("StopMuCluster");
    larcv::PSet stopmu_foxtrot_pset = tagger_pset.get<larcv::PSet>("StopMuFoxTrot");    
    larcv::PSet untagged_cluster_ps = tagger_pset.get<larcv::PSet>("ContainedGroupAlgo");
    larcv::PSet untagged_match_pset = tagger_pset.get<larcv::PSet>("TaggerFlashMatchAlgo");

    cfg.input_write_cfg  = tagger_pset.get<larcv::PSet>("InputWriteConfig");
    cfg.thrumu_write_cfg = tagger_pset.get<larcv::PSet>("ThruMuWriteConfig");
    cfg.stopmu_write_cfg = tagger_pset.get<larcv::PSet>("StopMuWriteConfig");
    cfg.croi_write_cfg   = tagger_pset.get<larcv::PSet>("CROIWriteConfig");

    cfg.sidetagger_cfg       = larlitecv::MakeBoundaryMuonTaggerAlgoConfigFromPSet( sidetagger_pset );
    cfg.flashtagger_cfg      = larlitecv::MakeFlashMuonTaggerAlgoConfigFromPSet( flashtagger_pset );
    cfg.thrumu_tracker_cfg   = larlitecv::ThruMuTrackerConfig::MakeFromPSet( thrumu_tracker_pset );
    cfg.stopmu_filterpts_cfg = larlitecv::MakeStopMuFilterSpacePointsConfigFromPSet( stopmu_filter_pset );
    cfg.stopmu_cluster_cfg   = larlitecv::makeStopMuClusterConfigFromPSet( stopmu_cluster_pset );
    cfg.stopmu_foxtrot_cfg   = larlitecv::StopMuFoxTrotConfig::makeFromPSet( stopmu_foxtrot_pset );
    cfg.untagged_cluster_cfg = larlitecv::ClusterGroupAlgoConfig::MakeClusterGroupAlgoConfigFromPSet( untagged_cluster_ps );
    cfg.croi_selection_cfg   = larlitecv::TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( untagged_match_pset );

    // run controls
    cfg.run_thrumu_tracker = tagger_pset.get<bool>("RunThruMuTracker",true);

    return cfg;
  }

}
