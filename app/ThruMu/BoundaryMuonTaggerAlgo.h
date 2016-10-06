#ifndef __LLCV_BOUNDARY_MUON_TAGGER_ALGO__
#define __LLCV_BOUNDARY_MUON_TAGGER_ALGO__

#include <vector>
#include <string>

#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"
#include "DataFormat/Pixel2DCluster.h"

#include "BoundaryMatchArrays.h"
#include "BoundaryEndPt.h"

namespace larlitecv {

  class ConfigBoundaryMuonTaggerAlgo {
    // container for boundary muon tagger algo configuration parameters.
    // fill them anyway you can
      

  public:
    ConfigBoundaryMuonTaggerAlgo() {};
    ~ConfigBoundaryMuonTaggerAlgo() {};

    bool checkOK() { return true; }; // dummy for now

  public:
    
    std::vector<float> thresholds;    ///< pixel threshold to count as a hit
    std::vector<int>   neighborhoods; ///< columns before and after to check for hits
    std::vector<int>   edge_win_wires; ///
    std::vector<int>   edge_win_times;
    std::vector<float> edge_win_hitthresh;
    std::vector<float> astar_thresholds; //< passed to astar config
  };

  class BoundaryMuonTaggerAlgo {
    // the algo

  public:

    BoundaryMuonTaggerAlgo() {};
    BoundaryMuonTaggerAlgo( ConfigBoundaryMuonTaggerAlgo& config ) {
      _config = config; //copy
    };
    virtual ~BoundaryMuonTaggerAlgo() {};

    enum { kOK=0, kErr_NotConfigured, kErr_BadInput };

  public:
    
    void configure( ConfigBoundaryMuonTaggerAlgo& cfg ) { _config = cfg; };
    void run() {};
    int searchforboundarypixels( const std::vector< larcv::Image2D >& imgs, // original image
				 std::vector< larcv::Image2D >& boundarypixelimgs ); // pixels consistent with boundary hits
    int clusterBoundaryPixels( const std::vector< larcv::Image2D >& imgs, // original image
			       const std::vector< larcv::Image2D >& matchedpixels, // pixels consistent with boundary hits
			       std::vector< std::vector<BoundaryEndPt> >& end_points ); // clustered end points on each plane
    int makePlaneTrackCluster( const larcv::Image2D& img, const larcv::Image2D& badchimg,
			       const std::vector< BoundaryEndPt >& top, const std::vector< BoundaryEndPt >& bot,
			       const std::vector< BoundaryEndPt >& upstream, const std::vector< BoundaryEndPt >& downstream,
			       const std::vector< BoundaryEndPt >& anode, const std::vector< BoundaryEndPt >& cathode,
			       std::vector< larcv::Pixel2DCluster >& trackclusters );

  protected:

    larlitecv::BoundaryMatchArrays m_matches;
    ConfigBoundaryMuonTaggerAlgo _config;

  };


};

#endif
