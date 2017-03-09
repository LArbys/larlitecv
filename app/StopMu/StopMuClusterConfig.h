#ifndef __STOPMUCLUSTERCONFIG__
#define __STOPMUCLUSTERCONFIG__

/*
  This class contains the adjustable parameters for the StopMuCluster class


 */

#include <vector>

// LArCV
#include "Base/PSet.h"

// larlitecv
#include "ThruMu/AStar3DAlgo.h"

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
    float link_stepsize;
    float astar_downsampling_factor;
    int num_passes;
    bool save_pass_images;
    bool dump_tagged_images;

    class PassConfig_t {
    public:
        float max_link_distance; //< maximum clusterl link distance
        float min_link_cosine; //< minimum cosine between links
        float max_extrema_row_diff; 
        float alldir_max_link_dist; //< if link is shorter than this dist, accept at all directions
        float max_extrema_triarea;
        AStar3DAlgoConfig astarcfg;
        PassConfig_t() {};
        virtual ~PassConfig_t() {};
    };
    std::vector< PassConfig_t > pass_configs;
    static PassConfig_t makePassConfigFromPSet( const larcv::PSet& pset );        
  };

  StopMuClusterConfig makeStopMuClusterConfigFromPSet( const larcv::PSet& pset );

}

#endif
