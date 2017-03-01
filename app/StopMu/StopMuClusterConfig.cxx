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
  }


  StopMuClusterConfig makeStopMuClusterConfigFromPSet( const larcv::PSet& pset ) {
    StopMuClusterConfig cfg;
    cfg.start_point_pixel_neighborhood = pset.get<int>("StartPointPixelNeighborhood");
    cfg.pixel_thresholds = pset.get< std::vector<float> >("PixelThresholds");
    cfg.dbscan_cluster_radius = pset.get< float >("ClusteringRadius");
    cfg.dbscan_cluster_minpoints = pset.get< int >("ClusteringMinPoints");
    cfg.max_link_distance = pset.get< float >("MaxLinkDistance");
    return cfg;
  }

}
