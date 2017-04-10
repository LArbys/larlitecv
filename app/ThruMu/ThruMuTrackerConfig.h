#ifndef __THRUMU_TRACKER_CONFIG__
#define __THRUMU_TRACKER_CONFIG__

/*
  ThruMuTrackerConfig

  Holds configurations for each of the passes made by the Tracker.

 */

// std lib
#include <vector>

// larcv
#include "Base/PSet.h"

// larlitecv
#include "Linear3DChargeTagger.h"

namespace larlitecv {

  class ThruMuTrackerConfig {
  public:
    ThruMuTrackerConfig(); //< loads with some defaults. not recomended.
    virtual ~ThruMuTrackerConfig() {};

    void setDefaults();

    static ThruMuTrackerConfig MakeFromPSet( const larcv::PSet& pset ); //< loads from pset. recommended.

    struct ThruMuPassConfig {
      bool run_linear_tagger; //< do we run the linearcharge3d algorithm this pass
      bool run_astar_tagger;  //< do we run the astar algorithm this pass
      int   linear3d_min_tracksize; //< how many steps should the line track have
      float linear3d_min_goodfraction; //< how many of the steps should see charge in all planes
      float linear3d_min_majoritychargefraction; //< how many should see charge in 2/3 planes
      float astar3d_min_goodfrac; //< what is the goodfrac the linear tagger should see before running astar
      float astar3d_min_majfrac;  //< what is the majfrac the linear tagger should see before running astar
      Linear3DChargeTaggerConfig linear3d_cfg; //< configuration of the linear 3d tagger
      AStar3DAlgoConfig astar3d_cfg; //< configuration of the astar algo config here
    };

    // Parameters
    
    int verbosity;
    int num_passes; //< number of passes
    std::vector< ThruMuPassConfig > pass_configs; //< configuration for each pass

    

  };

}

#endif
