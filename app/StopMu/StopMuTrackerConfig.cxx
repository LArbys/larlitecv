#include "StopMuTrackerConfig.h"

namespace larlitecv {

  StopMuTrackerConfig::StopMuTrackerConfig() {
    setDefaults();
  }

  void StopMuTrackerConfig::setDefaults() {
    // set some defaults
    nplanes = 3;
    skeleton_kernel_size = 3;
    start_point_pixel_neighborhood = 10;
    pixel_thresholds.resize(nplanes,10.0);
    dbscan_cluster_radius = 50.0;
    dbscan_cluster_minpoints = 5;
    trigger_tick = 3200;
    usec_per_tick = 0.5;
    cm_per_wire = 0.3;
    fitter_distance_weight = 2.0;
    fitter_bend_weight = 0.1;
    fitter_charge_weight = 0.;
    fitter_step_size_cm  = 2.0;
    max_steppoint_dist_from_cluster = 0.3*20.0;
  }


}
