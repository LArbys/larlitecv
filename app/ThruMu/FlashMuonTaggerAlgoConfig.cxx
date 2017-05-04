#include "FlashMuonTaggerAlgoConfig.h"

namespace larlitecv {

  void FlashMuonTaggerAlgoConfig::setdefaults() {

    pixel_value_threshold.resize(3,10.0);
    clustering_minpoints.resize(3,3);
    clustering_radius.resize(3,5.0);
    endpoint_time_neighborhood.resize(3,10);
    verbosity = 0;
    trigger_tick = 3200;
    usec_per_tick = 0.5;
    drift_distance = 250.0;
    drift_velocity = 0.114;
    flash_zrange_extension = 1.1;
    max_triarea = 10.0;
    max_nsegments_per_flash = 3;
    cathode_drift_tick_correction = 240.0;
    endpoint_clustering_algo = "segments";
    
  }


  FlashMuonTaggerAlgoConfig MakeFlashMuonTaggerAlgoConfigFromPSet( const larcv::PSet& flashtagger_pset ) {
    FlashMuonTaggerAlgoConfig flashtagger_cfg;
    flashtagger_cfg.pixel_value_threshold        = flashtagger_pset.get< std::vector<float> >( "ChargeThreshold" ); // one for each plane
    flashtagger_cfg.clustering_minpoints         = flashtagger_pset.get< std::vector<int> >( "ClusteringMinPoints" );
    flashtagger_cfg.clustering_radius            = flashtagger_pset.get< std::vector<double> >( "ClusteringRadius" );
    flashtagger_cfg.endpoint_time_neighborhood   = flashtagger_pset.get< std::vector<int> >( "EndpointTimeNeighborhood" );
    flashtagger_cfg.verbosity                    = flashtagger_pset.get< int >( "Verbosity", 2 );
    flashtagger_cfg.trigger_tick                 = flashtagger_pset.get< float >( "TriggerTick", 3200.0 );
    flashtagger_cfg.usec_per_tick                = flashtagger_pset.get< float >( "MicrosecondsPerTick", 0.5 );
    flashtagger_cfg.drift_distance               = flashtagger_pset.get< float >( "DriftDistance", 250.0 );
    flashtagger_cfg.drift_velocity               = flashtagger_pset.get< float >( "DriftVelocity", 0.114 );
    flashtagger_cfg.flash_zrange_extension       = flashtagger_pset.get< float >( "FlashZRangeExtension", 1.1 );
    flashtagger_cfg.max_triarea                  = flashtagger_pset.get< float >( "MaxTriArea", 10.0 );
    flashtagger_cfg.max_nsegments_per_flash      = flashtagger_pset.get< int >( "MaxNumSegmentsPerFlash", 3 );
    flashtagger_cfg.cathode_drift_tick_correction = flashtagger_pset.get< float >( "CathodeDriftTickCorrection" );
    flashtagger_cfg.endpoint_clustering_algo     = flashtagger_pset.get< std::string >("EndPointClusteringAlgo");
    return flashtagger_cfg;
  }

  
}
