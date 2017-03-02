#ifndef __STOPMUCLUSTERCONFIG__
#define __STOPMUCLUSTERCONFIG__

/*
  This class contains the adjustable parameters for the StopMuCluster class


 */

#include <vector>

#include "Base/PSet.h"

namespace larlitecv {

  class StopMuClusterConfig {

  public:
    StopMuClusterConfig();
    virtual ~StopMuClusterConfig() {};

    void setDefaults();
    
    int verbosity; //< verbosity. 0=most quiet

    int start_point_pixel_neighborhood; //< pixels around candidate starting points that cannot be vetoed by thrumu tagger
    std::vector<float> pixel_thresholds; //< threshold of charge to consider pixel interesting or not
    float dbscan_cluster_radius; //< pixels have to be within this radius to be clustered together
    int dbscan_cluster_minpoints; //< minimum number of points in a cluster
    float max_link_distance; //< maximum clusterl link distance
    float min_link_cosine; //< minimum cosine between links
    float max_extrema_row_diff; 
    float max_extrema_triarea;
    float link_stepsize;
    float astar_downsampling_factor;

  };

  StopMuClusterConfig makeStopMuClusterConfigFromPSet( const larcv::PSet& pset );

}

#endif
