#include "FlashMuonTaggerConfig.h"

namespace larlitecv {

  void FlashMuonTaggerConfig::setdefaults() {

    pixel_value_threshold.resize(3,10.0);
    search_row_radius = 10;
    clustering_time_neighborhood.resize(3,20);
    clustering_wire_neighborhood.resize(3,20);
    clustering_minpoints.resize(3,3);
    clustering_radius.resize(3,2.0);
    endpoint_time_neighborhood.resize(3,10);
    verbosity = 2;
    trigger_tick = 3200;
    usec_per_tick = 0.5;
    drift_distance = 250.0;
    drift_velocity = 0.114;
    flash_zrange_extension = 1.1;
    max_triarea = 10.0;
    max_triarea_tight  = 1.0;
    
  }


  FlashMuonTaggerConfig MakeFlashMuonTaggerConfigFromPSet( const larcv::PSet& flashtagger_pset ) {
    FlashMuonTaggerConfig flashtagger_cfg;
    
    flashtagger_cfg.pixel_value_threshold        = flashtagger_pset.get< std::vector<float> >( "ChargeThreshold" ); // one for each plane
    flashtagger_cfg.search_row_radius            = flashtagger_pset.get< int >( "SearchRowRadius" ); // one for each plane
    flashtagger_cfg.clustering_time_neighborhood = flashtagger_pset.get< std::vector<int> >( "ClusteringTimeWindow" );
    flashtagger_cfg.clustering_wire_neighborhood = flashtagger_pset.get< std::vector<int> >( "ClusteringWireWindow" );
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
    flashtagger_cfg.max_triarea_tight            = flashtagger_pset.get< float >( "MaxTriAreaTight", 1.0 );
    
    return flashtagger_cfg;
  }

  
}
