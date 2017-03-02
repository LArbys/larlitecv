#ifndef __DBCLUSTER_MERGER__
#define __DBCLUSTER_MERGER__

// LArCV
#include "dbscan/DBSCANAlgo.h"
#include "DataFormat/Image2D.h"

namespace larlitecv {

  // parameters to control this algorithm
  struct DBClusterMergerPars_t {
    int astar_image_padding;     // sets the extra rows and columns around the clusters of interest to consider
    int astar_neighborhood;      // sets size of neighborhood that A* algorithm considers when connecting one pixel to another
    int astar_start_pad;         // region around start points that A* considers OK, even if no charge. allows sloppy start point
    int astar_end_pad;           // region around start points that A* considers OK, even if no charge. allows sloppy end point
    float astar_pixel_threshold; // value pixel (on non-badch) must have if path allowed to continue through it
    DBClusterMergerPars_t() {
      astar_pixel_threshold = 10.0; // one should really set this as it depends on compression factor of image			
      // defaults that should be workable
      astar_image_padding = 10;
      astar_neighborhood = 1;
      astar_start_pad = 2;
      astar_end_pad = 2;
    };
  };

  class DBClusterMergerAlgo {
  public:

    DBClusterMergerAlgo() {};
    virtual ~DBClusterMergerAlgo() {};

    static bool willClustersMergeThroughBadchs( const int clusterid_a, const int clusterid_b, 
      const dbscan::dbPoints& pixels, const dbscan::dbscanOutput& cluster_info, 
      const larcv::Image2D& img, const larcv::Image2D& badch_img, const DBClusterMergerPars_t& config );

  };

}
#endif