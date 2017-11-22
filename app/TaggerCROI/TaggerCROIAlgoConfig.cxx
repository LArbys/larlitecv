#include "TaggerCROIAlgoConfig.h"

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite
#include "FhiclLite/FhiclLiteUtilFunc.h"

namespace larlitecv {

  TaggerCROIAlgoConfig::TaggerCROIAlgoConfig()
    : input_write_cfg(larcv::PSet("input")),
      thrumu_write_cfg(larcv::PSet("thrumu")),
      stopmu_write_cfg(larcv::PSet("stopmu")),
      croi_write_cfg(larcv::PSet("croi"))
  {
    run_thrumu_tracker = true;
    recluster_stop_and_thru = true;
    use_truth_endpoints = false;
    use_truth_muontracks = false;    
    verbosity = 0;
  }

  TaggerCROIAlgoConfig TaggerCROIAlgoConfig::makeConfigFromFile( std::string filepath ) {

    TaggerCROIAlgoConfig cfg;

    larcv::PSet root_pset                = larcv::CreatePSetFromFile( filepath );
    larcv::PSet tagger_pset              = root_pset.get<larcv::PSet>("TaggerCROI");
    larcv::PSet sidetagger_pset          = tagger_pset.get<larcv::PSet>("BMTSideTagger");
    larcv::PSet flashtagger_pset         = tagger_pset.get<larcv::PSet>("BMTFlashTagger");
    larcv::PSet thrumu_tracker_pset      = tagger_pset.get<larcv::PSet>("ThruMuTracker");
    larcv::PSet stopmu_filter_pset       = tagger_pset.get<larcv::PSet>("StopMuSpacePointsFilter");
    larcv::PSet stopmu_cluster_pset      = tagger_pset.get<larcv::PSet>("StopMuCluster");
    larcv::PSet stopmu_foxtrot_pset      = tagger_pset.get<larcv::PSet>("StopMuFoxTrot");
    larcv::PSet untagged_cluster_ps      = tagger_pset.get<larcv::PSet>("ContainedGroupAlgo");
    larcv::PSet untagged_match_pset      = tagger_pset.get<larcv::PSet>("TaggerFlashMatchAlgo");
    //larcv::PSet general_flash_match_pset = tagger_pset.get<larcv::PSet>("GeneralFlashMatchAlgo"); // new for the general flash matching occurring between the tracks and the flashes of light.

    cfg.input_write_cfg  = tagger_pset.get<larcv::PSet>("InputWriteConfig");
    cfg.thrumu_write_cfg = tagger_pset.get<larcv::PSet>("ThruMuWriteConfig");
    cfg.stopmu_write_cfg = tagger_pset.get<larcv::PSet>("StopMuWriteConfig");
    cfg.croi_write_cfg   = tagger_pset.get<larcv::PSet>("CROIWriteConfig");

    cfg.sidetagger_cfg          = larlitecv::MakeBoundaryMuonTaggerAlgoConfigFromPSet( sidetagger_pset );
    cfg.flashtagger_cfg         = larlitecv::MakeFlashMuonTaggerAlgoConfigFromPSet( flashtagger_pset );
    cfg.thrumu_tracker_cfg      = larlitecv::ThruMuTrackerConfig::MakeFromPSet( thrumu_tracker_pset );
    cfg.stopmu_filterpts_cfg    = larlitecv::MakeStopMuFilterSpacePointsConfigFromPSet( stopmu_filter_pset );
    cfg.stopmu_cluster_cfg      = larlitecv::makeStopMuClusterConfigFromPSet( stopmu_cluster_pset );
    cfg.stopmu_foxtrot_cfg      = larlitecv::StopMuFoxTrotConfig::makeFromPSet( stopmu_foxtrot_pset );
    cfg.untagged_cluster_cfg    = larlitecv::ClusterGroupAlgoConfig::MakeClusterGroupAlgoConfigFromPSet( untagged_cluster_ps );
    cfg.croi_selection_cfg      = larlitecv::TaggerFlashMatchAlgoConfig::FromPSet( untagged_match_pset );
    //cfg.general_flash_match_cfg = larlitecv::GeneralFlashMatchAlgoConfig::MakeConfigFromPSet( general_flash_match_pset );

    // run controls
    cfg.run_thrumu_tracker      = tagger_pset.get<bool>("RunThruMuTracker",true);
    cfg.recluster_stop_and_thru = tagger_pset.get<bool>("ReclusterStopAndThruMu",true);
    cfg.use_truth_endpoints     = tagger_pset.get<bool>("UseTruthEndPoints",false);
    cfg.use_truth_muontracks    = tagger_pset.get<bool>("UseTruthMuonTracks",false);    
    cfg.verbosity               = tagger_pset.get<int>("Verbosity");

    cfg.larcv_image_producer    = tagger_pset.get<std::string>("LArCVImageProducer");
    cfg.larcv_chstatus_producer = tagger_pset.get<std::string>("LArCVChStatusProducer");
    cfg.DeJebWires              = tagger_pset.get<bool>("DeJebWires");
    cfg.jebwiresfactor          = tagger_pset.get<float>("JebWiresFactor");
    cfg.emptych_thresh          = tagger_pset.get< std::vector<float> >("EmptyChannelThreshold");
    cfg.chstatus_datatype       = tagger_pset.get<std::string>("ChStatusDataType");
    cfg.opflash_producers       = tagger_pset.get< std::vector<std::string> >( "OpflashProducers" );
    cfg.trigger_producer        = tagger_pset.get< std::string >( "TriggerProducer" );    
    cfg.RunThruMu               = tagger_pset.get<bool>("RunThruMu");
    cfg.RunStopMu               = tagger_pset.get<bool>("RunStopMu");
    cfg.RunCROI                 = tagger_pset.get<bool>("RunCROI");
    cfg.save_thrumu_space       = tagger_pset.get<bool>("SaveThruMuSpace", false);
    cfg.save_stopmu_space       = tagger_pset.get<bool>("SaveStopMuSpace", false);
    cfg.save_croi_space         = tagger_pset.get<bool>("SaveCROISpace",   false);
    cfg.save_mc                 = tagger_pset.get<bool>("SaveMC",false);
    cfg.load_mctrack            = tagger_pset.get<bool>("LoadMCTrack",false);
    cfg.mctrack_producer        = tagger_pset.get<std::string>("MCTrackProducer","mcreco");        
    cfg.skip_empty_events       = tagger_pset.get<bool>("SkipEmptyEvents",false);
    cfg.apply_unipolar_hack     = tagger_pset.get<bool>("ApplyUnipolarHack",false);

    return cfg;
  }

}
