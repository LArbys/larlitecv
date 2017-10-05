#include "StopMuClusterConfig.h"

namespace larlitecv {


  StopMuClusterConfig::StopMuClusterConfig() {
    setDefaults();
  }

  void StopMuClusterConfig::setDefaults() {
    // set some defaults
    start_point_pixel_neighborhood = 5;
    pixel_thresholds.resize(3,10.0);
    dbscan_cluster_radius = 10.0;
    dbscan_cluster_minpoints = 10;
    link_stepsize = 0.3;
    astar_downsampling_factor = 4; 
    num_passes = 2;
    save_pass_images = false;
    dump_tagged_images = false;
    PassConfig_t passcfg;    
    passcfg.max_link_distance = 150;
    passcfg.min_link_cosine = 0.8;
    passcfg.alldir_max_link_dist = 20.0;
    passcfg.max_extrema_row_diff = 2.0;
    passcfg.max_extrema_triarea = 2.0;
    passcfg.astarcfg.astar_threshold.resize(3,10.0);
    passcfg.astarcfg.astar_neighborhood.resize(3,4);
    passcfg.astarcfg.astar_start_padding = 4;
    passcfg.astarcfg.astar_end_padding = 4;    
    passcfg.astarcfg.lattice_padding = 10;
    passcfg.astarcfg.accept_badch_nodes = true;
    passcfg.astarcfg.min_nplanes_w_hitpixel = 3;
    passcfg.astarcfg.restrict_path = false;
    passcfg.astarcfg.path_restriction_radius = 10.0;
    passcfg.astarcfg.verbosity = 0;

    // setup up two identical passes
    pass_configs.push_back( passcfg );
    pass_configs.push_back( passcfg );    

  }

  StopMuClusterConfig::PassConfig_t StopMuClusterConfig::makePassConfigFromPSet( const larcv::PSet& pass_pset ) {
    PassConfig_t passcfg;
    passcfg.max_link_distance    = pass_pset.get< float >("MaxLinkDistance");
    passcfg.min_link_cosine      = pass_pset.get< float >("MinLinkCosine");
    passcfg.alldir_max_link_dist = pass_pset.get< float >("AllDirMaxLinkDistance");
    passcfg.max_extrema_row_diff = pass_pset.get< float >("MaxExtremaRowDiff");
    passcfg.max_extrema_triarea  = pass_pset.get< float >("MaxExtremaTriArea");
    passcfg.astarcfg             = larcv::AStar3DAlgoConfig::MakeFromPSet( pass_pset.get<larcv::PSet>("AStarConfig") );
    return passcfg;
  }

  StopMuClusterConfig makeStopMuClusterConfigFromPSet( const larcv::PSet& pset ) {

    StopMuClusterConfig cfg;
    cfg.start_point_pixel_neighborhood = pset.get<int>("StartPointPixelNeighborhood");
    cfg.pixel_thresholds = pset.get< std::vector<float> >("PixelThresholds");
    cfg.num_passes = pset.get<int>("NumPasses");
    cfg.dbscan_cluster_radius = pset.get< float >("ClusteringRadius");
    cfg.dbscan_cluster_minpoints = pset.get< int >("ClusteringMinPoints");
    cfg.link_stepsize = pset.get<float>("LinkStepSize",0.3);
    cfg.astar_downsampling_factor = pset.get<float>("AStarDownsamplingFactor",4.0);
    cfg.save_pass_images = pset.get<bool>("SavePassImages");
    cfg.dump_tagged_images = pset.get<bool>("DumpTaggedImages");    

    cfg.pass_configs.clear();
    for (int ipass=0; ipass<cfg.num_passes; ipass++) {
      char zpass[20];
      sprintf(zpass, "Pass%d",ipass+1);
      larcv::PSet pass_pset = pset.get<larcv::PSet>(std::string(zpass));
      StopMuClusterConfig::PassConfig_t passcfg = StopMuClusterConfig::makePassConfigFromPSet(pass_pset);        
      cfg.pass_configs.emplace_back( std::move(passcfg) );
    }

    return cfg;
  }

}
