#ifndef __STOPMUTRACKERCONFIG__
#define __STOPMUTRACKERCONFIG__

/*
  This class contains the adjustable parameters for the StopMuTracker class


 */

#include <vector>

namespace larlitecv {

  class StopMuTrackerConfig {

  public:
    StopMuTrackerConfig();
    virtual ~StopMuTrackerConfig() {};

    void setDefaults();
    void checkParameters() {};
    
    int verbosity; //< verbosity. 0=most quiet

    int nplanes; //< number of planes
    int skeleton_kernel_size; //< kernel size for skeletonization operator
    int start_point_pixel_neighborhood; //< pixels around candidate starting points that cannot be vetoed by thrumu tagger
    std::vector<float> pixel_thresholds; //< threshold of charge to consider pixel interesting or not
    float dbscan_cluster_radius; //< pixels have to be within this radius to be clustered together
    int dbscan_cluster_minpoints; //< minimum number of points in a cluster
    int trigger_tick; // tick number where trigger is expected to occur
    float usec_per_tick; // microseconds per TPC tick
    float cm_per_wire;   // distance between wires
    float fitter_bend_weight;
    float fitter_distance_weight;
    float fitter_charge_weight;
    float fitter_step_size_cm;
    float max_steppoint_dist_from_cluster;
  };
}

#endif
