#include "ConfigBoundaryMuonTaggerAlgo.h"

namespace larlitecv {

  void ConfigBoundaryMuonTaggerAlgo::setdefaults() {

    save_endpt_images = false; 
    hitsearch_uses_badchs = false;
    ticks_per_full_drift = 4650;
    type_modifier.resize(4,1.0);
    type_modifier[2] = 0.5;
    type_modifier[3] = 0.5;

    neighborhoods.resize(3,10);
    thresholds.resize(3,10.0);
    emptych_thresh.resize(10.0);
    edge_win_wires.resize(12,20);
    edge_win_times.resize(12,20);
    edge_win_hitthresh.resize(12,10.0);
    boundary_cluster_minpixels.resize(3,10);
    boundary_cluster_radius.resize(3,10);
    astar_thresholds.resize(3,10.0);
    astar_neighborhood.resize(3,10.0);
    tag_neighborhood.resize(3,2);
    verbosity = 0;
    linear3d_min_tracksize = 5;
    linear3d_min_goodfraction = 0.95;
    linear3d_min_majoritychargefraction = 0.5;
    astar3d_min_goodfrac = 0.2;
    astar3d_min_majfrac = 0.2;
  }

  ConfigBoundaryMuonTaggerAlgo MakeConfigBoundaryMuonTaggerAlgoFromPSet( const larcv::PSet& pset ) {

    ConfigBoundaryMuonTaggerAlgo config;

    config.neighborhoods              = pset.get< std::vector<int> >("Neighborhoods");
    config.thresholds                 = pset.get< std::vector<float> >( "Thresholds" );
    config.emptych_thresh             = pset.get< std::vector<float> >( "EmptyChannelThrehsold" );
    config.edge_win_wires             = pset.get< std::vector<int> >( "EdgeWinWires" );
    config.edge_win_times             = pset.get< std::vector<int> >( "EdgeWinTimes" );
    config.edge_win_hitthresh         = pset.get< std::vector<float> >( "EdgeWinHitThreshold" );
    config.boundary_cluster_minpixels = pset.get< std::vector<int> >( "BoundaryClusterMinPixels" );
    config.boundary_cluster_radius    = pset.get< std::vector<float> >( "BoundaryClusterRadius" );
    config.astar_thresholds           = pset.get< std::vector<float> >( "AStarThresholds" );
    config.astar_neighborhood         = pset.get< std::vector<int> >( "AStarNeighborhood" );
    config.tag_neighborhood           = pset.get< std::vector<int> >("TaggingNeighborhood");    
    config.save_endpt_images          = pset.get<bool>("SaveMatchImages",false);
    config.hitsearch_uses_badchs      = pset.get<bool>("UseBadChannels",true);
    config.ticks_per_full_drift       = pset.get<float>("TicksPerFullDrift",4650.0);
    config.verbosity                  = pset.get<int>("Verbosity",0);
    config.astar_cfg                  = larlitecv::AStar3DAlgoConfig::MakeFromPSet( pset.get< larcv::PSet >( "AStarConfig" ) );
    config.linear3d_cfg               = larlitecv::Linear3DFitterConfig::makeFromPSet( pset.get< larcv::PSet >( "Linear3DConfig" ) );
    config.linear3d_min_tracksize     = pset.get<int>("Linear3DMinTrackSize");
    config.linear3d_min_goodfraction  = pset.get<float>("Linear3DMinGoodFraction");
    config.linear3d_min_majoritychargefraction = pset.get<float>("Linear3DMinMajorityChargeFraction");
    config.astar3d_min_goodfrac       = pset.get<float>("AStar3DMinGoodFraction");
    config.astar3d_min_majfrac        = pset.get<float>("AStar3DMinMajorityChargeFraction");
               
    return config;
  }

}
