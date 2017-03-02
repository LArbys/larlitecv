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
    dbscan_cluster_minpoints = 30;
    max_link_distance = 150;
    min_link_cosine = 0.8;
    max_extrema_row_diff = 2.0;
    max_extrema_triarea = 2.0;
    link_stepsize = 0.3;
    astar_downsampling_factor = 4;
  }


  StopMuClusterConfig makeStopMuClusterConfigFromPSet( const larcv::PSet& pset ) {
    StopMuClusterConfig cfg;
    cfg.start_point_pixel_neighborhood = pset.get<int>("StartPointPixelNeighborhood");
    cfg.pixel_thresholds = pset.get< std::vector<float> >("PixelThresholds");
    cfg.dbscan_cluster_radius = pset.get< float >("ClusteringRadius");
    cfg.dbscan_cluster_minpoints = pset.get< int >("ClusteringMinPoints");
    cfg.max_link_distance = pset.get< float >("MaxLinkDistance");
    cfg.min_link_cosine = pset.get< float >("MinLinkCosine");
    cfg.max_extrema_row_diff = pset.get< float >("MaxExtremaRowDiff");
    cfg.max_extrema_triarea = pset.get< float >("MaxExtremaTriArea");
    return cfg;
  }

}
