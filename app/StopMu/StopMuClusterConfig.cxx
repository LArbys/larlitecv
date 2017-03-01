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
  }


}
