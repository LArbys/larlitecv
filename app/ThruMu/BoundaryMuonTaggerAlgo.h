#ifndef __LLCV_BOUNDARY_MUON_TAGGER_ALGO__
#define __LLCV_BOUNDARY_MUON_TAGGER_ALGO__

#include <vector>
#include <string>

#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"
#include "DataFormat/Pixel2DCluster.h"

#include "PMTWeights/WireData.h"

#include "BoundaryMatchArrays.h"
#include "BoundaryEndPt.h"
#include "BMTrackCluster2D.h"
#include "BMTrackCluster3D.h"

namespace larlitecv {

  class ConfigBoundaryMuonTaggerAlgo {
    // container for boundary muon tagger algo configuration parameters.
    // fill them anyway you can
      

  public:
    ConfigBoundaryMuonTaggerAlgo() {};
    ~ConfigBoundaryMuonTaggerAlgo() {};

    bool checkOK() { return true; }; // dummy for now

  public:

    std::vector<float> emptych_thresh; ///< pixel thresholds below which if wire stays, it is marked as empty (value per plane)
    std::vector<float> thresholds;    ///< pixel threshold to count as a hit
    std::vector<int>   neighborhoods; ///< columns before and after to check for hits
    std::vector<int>   edge_win_wires; ///
    std::vector<int>   edge_win_times;
    std::vector<float> edge_win_hitthresh;
    std::vector<int>   boundary_cluster_minpixels; ///< min number of pixels for dbscan clustering of boundary pixels
    std::vector<float>   boundary_cluster_radius;    ///< nearest neighbor radius for dbscan clustering of boundary pixels
    std::vector<float> astar_thresholds; //< passed to astar config
    std::vector<int>   astar_neighborhood; //< passed to astar config
  };

  class BoundaryMuonTaggerAlgo {
    // the algo

  public:

    BoundaryMuonTaggerAlgo() { loadGeoInfo(); };
    BoundaryMuonTaggerAlgo( ConfigBoundaryMuonTaggerAlgo& config ) {
      _config = config; //copy
      loadGeoInfo();
    };
    virtual ~BoundaryMuonTaggerAlgo() {};

    enum { kOK=0, kErr_NotConfigured, kErr_BadInput };

  public:
    
    void configure( ConfigBoundaryMuonTaggerAlgo& cfg ) { _config = cfg; };
    void run();
    int searchforboundarypixels( const std::vector< larcv::Image2D >& imgs, // original image
				 const std::vector< larcv::Image2D >& badchs, // image with bad channels marked
				 std::vector< larcv::Image2D >& boundarypixelimgs ); // pixels consistent with boundary hits
    int clusterBoundaryPixels( const std::vector< larcv::Image2D >& imgs, // original image
			       const std::vector< larcv::Image2D >& matchedpixels, // pixels consistent with boundary hits
			       std::vector< std::vector<BoundaryEndPt> >& end_points ); // clustered end points on each plane
    int searchforboundarypixels3D( const std::vector< larcv::Image2D >& imgs, // original image
				   const std::vector< larcv::Image2D >& badchs, // image with bad channels marked
				   std::vector< larcv::Image2D >& boundarypixelimgs ); // pixels consistent with boundary hits
    int clusterBoundaryPixels3D( const std::vector< larcv::Image2D >& matchedpixels, // pixels consistent with boundary hits
				 std::vector< std::vector<BoundaryEndPt> >& end_points ); // list of end point triples
    int makePlaneTrackCluster( const larcv::Image2D& img, const larcv::Image2D& badchimg,
			       const std::vector< BoundaryEndPt >& top, const std::vector< BoundaryEndPt >& bot,
			       const std::vector< BoundaryEndPt >& upstream, const std::vector< BoundaryEndPt >& downstream,
			       const std::vector< BoundaryEndPt >& anode, const std::vector< BoundaryEndPt >& cathode,
			       const std::vector< BoundaryEndPt >& imgends,
			       std::vector< larlitecv::BMTrackCluster2D >& trackclusters );
    int markImageWithTrackClusters( const std::vector<larcv::Image2D>& imgs, const std::vector< std::vector< larlitecv::BMTrackCluster2D > >& trackclusters,
				    std::vector<larcv::Image2D>& markedimgs );
    void matchTracksStage1( const std::vector< larcv::Image2D >& imgs, const std::vector< std::vector< larlitecv::BMTrackCluster2D >* >& plane2dtracks, 
			    std::vector< larlitecv::BMTrackCluster3D >& output  );
    bool doTracksMatch( const larlitecv::BMTrackCluster2D& track1, const larlitecv::BMTrackCluster2D& track2, 
			float& start_t_diff, float& end_t_diff, bool& start2start );
			

    // Wire Geometry info
    void loadGeoInfo();
    std::map<int,larcv::pmtweights::WireData> m_WireData; // key is plane ID, value is class with wire info    
    int fNPMTs;
    float pmtpos[32][3];

  protected:


    larlitecv::BoundaryMatchArrays m_matches;
    ConfigBoundaryMuonTaggerAlgo _config;

  };


};

#endif
