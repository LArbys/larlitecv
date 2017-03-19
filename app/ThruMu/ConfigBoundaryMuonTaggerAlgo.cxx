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
    track_step_size = 0.3;
  }

  void ConfigBoundaryMuonTaggerAlgo::print() {
    std::cout << "==================================================" << std::endl;
    std::cout << "ConfigBoundaryMuonTaggerAlgo" << std::endl;

    std::cout << " Pixel Thresholds: "; for ( size_t i=0; i<thresholds.size(); i++ ) std::cout << thresholds[i] << " "; std::cout << std::endl;
    std::cout << " Search Neighborhood: "; for ( size_t i=0; i<neighborhoods.size(); i++ ) std::cout << neighborhoods[i] << " "; std::cout << std::endl; 
    std::cout << " Hit Search Uses BadChs: " << hitsearch_uses_badchs << std::endl;
    std::cout << " Boundary Cluster Min Pixels: "; for ( size_t i=0; i<boundary_cluster_minpixels.size(); i++ ) std::cout << boundary_cluster_minpixels[i] << " "; std::cout << std::endl;
    std::cout << " Boundary Cluster Radius: "; for ( size_t i=0; i<boundary_cluster_radius.size(); i++ ) std::cout << boundary_cluster_radius[i] << " "; std::cout << std::endl;
    std::cout << " Tag Neighborhood: "; for ( size_t i=0; i<tag_neighborhood.size(); i++ ) std::cout << tag_neighborhood[i] << " "; std::cout << std::endl;
    std::cout << " Type Modifer: "; for ( size_t i=0; i<type_modifier.size(); i++ ) std::cout << type_modifier[i] << " "; std::cout << std::endl;
    std::cout << " Linear3D Min Majority Charge Fraction: " << linear3d_min_majoritychargefraction << std::endl;
    std::cout << " Linear3D Min Good Fraction: " << linear3d_min_goodfraction << std::endl;
    std::cout << " Linear3D Min Track Size: " << linear3d_min_tracksize << std::endl;
    std::cout << " AStar3D Min Good Fraction: " << astar3d_min_goodfrac << std::endl;
    std::cout << " AStar3D Min Majority Fraction: " << astar3d_min_majfrac << std::endl;
    std::cout << " TicksPerFullDrift: " << ticks_per_full_drift << std::endl;
    std::cout << " Verbosity: " << verbosity << std::endl;
    std::cout << " Save Endpoint Images: " << save_endpt_images << std::endl;
    std::cout << "==================================================" << std::endl;    
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

    config.tag_neighborhood           = pset.get< std::vector<int> >("TaggingNeighborhood");    
    config.save_endpt_images          = pset.get<bool>("SaveMatchImages",false);
    config.hitsearch_uses_badchs      = pset.get<bool>("UseBadChannels",true);
    config.ticks_per_full_drift       = pset.get<float>("TicksPerFullDrift",4650.0);
    config.verbosity                  = pset.get<int>("Verbosity",0);
    config.astar_cfg                  = larlitecv::AStar3DAlgoConfig::MakeFromPSet( pset.get< larcv::PSet >( "AStarConfig" ) );
    config.linear3d_cfg               = larlitecv::Linear3DChargeTaggerConfig::makeFromPSet( pset.get< larcv::PSet >( "Linear3DConfig" ) );
    config.linear3d_min_tracksize     = pset.get<int>("Linear3DMinTrackSize");
    config.linear3d_min_goodfraction  = pset.get<float>("Linear3DMinGoodFraction");
    config.linear3d_min_majoritychargefraction = pset.get<float>("Linear3DMinMajorityChargeFraction");
    config.astar3d_min_goodfrac       = pset.get<float>("AStar3DMinGoodFraction");
    config.astar3d_min_majfrac        = pset.get<float>("AStar3DMinMajorityChargeFraction");
    config.track_step_size            = pset.get<float>("TrackStepSize");
               
    return config;
  }

}
