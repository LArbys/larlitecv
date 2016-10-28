#ifndef __LLCV_BOUNDARY_MUON_TAGGER_ALGO__
#define __LLCV_BOUNDARY_MUON_TAGGER_ALGO__

#include <vector>
#include <string>

// larcv
#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"
#include "DataFormat/Pixel2DCluster.h"
#include "dbscan/DBSCANAlgo.h"
#include "UBWireTool/WireData.h"

#include "BoundaryMatchArrays.h"
#include "BoundaryMatchAlgo.h"
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
    virtual ~BoundaryMuonTaggerAlgo();

    enum { kOK=0, kErr_NotConfigured, kErr_BadInput };

  public:
    
    void configure( ConfigBoundaryMuonTaggerAlgo& cfg ) { _config = cfg; };
    void run();
    int searchforboundarypixels3D( const std::vector< larcv::Image2D >& imgs, // original image
				   const std::vector< larcv::Image2D >& badchs, // image with bad channels marked
				   std::vector< std::vector<BoundaryEndPt> >& end_points, ///list of end point triples
				   std::vector< larcv::Image2D >& boundarypixelimgs, // pixels consistent with boundary hits
				   std::vector< larcv::Image2D >& boundaryspaceptsimgs ); // points in real-space consistent with boundary hits
    int makeTrackClusters3D( std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
			     const std::vector< const std::vector< BoundaryEndPt >* >& spacepts,
			     std::vector< std::vector< larlitecv::BMTrackCluster2D > >& trackclusters );
    int markImageWithTrackClusters( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
				    const std::vector< std::vector< larlitecv::BMTrackCluster2D > >& trackclusters, std::vector<int>& goodlist,
				    std::vector<larcv::Image2D>& markedimgs );
    BMTrackCluster2D runAstar( const BoundaryEndPt& start, const BoundaryEndPt& end, const larcv::Image2D& img, const larcv::Image2D& badchimg,
			       int start_pad, int end_pad, int verbose=2, bool use_badchs=false );
    void process2Dtracks( std::vector< std::vector< larlitecv::BMTrackCluster2D > >& trackclusters2D,
			  const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
			  std::vector< BMTrackCluster3D >& tracks, std::vector<int>& goodlist );

			

    // Wire Geometry info
    int fNPMTs;
    float pmtpos[32][3];

  protected:

    void loadGeoInfo();

    ConfigBoundaryMuonTaggerAlgo _config;

    BoundaryMatchAlgo* matchalgo_tight;
    BoundaryMatchAlgo* matchalgo_loose;
    
    void getClusterEdges( const dbscan::dbPoints& points, const std::vector< larcv::Image2D >& imgs, 
			  const dbscan::dbscanOutput& clout, int idx_cluster,
			  int& idxhit_tmin, int& idxhit_tmax, int& idxhit_wmin, int& idxhit_wmax );

    BMTrackCluster3D process2Dtrack( std::vector< larlitecv::BMTrackCluster2D >& track2d, 
				     const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v );
    bool compare2Dtrack( const std::vector< BMTrackCluster2D >& track2d, const BMTrackCluster3D& track3d, const larcv::ImageMeta& meta,
			 float path_radius_cm, float endpt_radius_cm );


  };


};

#endif
