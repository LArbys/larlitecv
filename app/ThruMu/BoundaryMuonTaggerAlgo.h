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

#include "ConfigBoundaryMuonTaggerAlgo.h"
#include "BoundaryMatchArrays.h"
#include "BoundaryMatchAlgo.h"
#include "BoundaryEndPt.h"
#include "BMTrackCluster2D.h"
#include "BMTrackCluster3D.h"

namespace larlitecv {

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
    typedef enum { top=0, bot, upstream, downstream, kNumCrossings } Crossings_t;

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

    class ClusterExtrema_t {
    public:
      std::vector<int> l; // left-most point
      std::vector<int> r; // right-most point
      std::vector<int> t; // top-most point
      std::vector<int> b; // bottom-most point
      int filled;
      ClusterExtrema_t() {
        l.resize(2,0);
        r.resize(2,0);
        t.resize(2,0);
        b.resize(2,0);
        filled = 0;
      };
      ~ClusterExtrema_t() {};
    };    

  protected:

    void loadGeoInfo();

    ConfigBoundaryMuonTaggerAlgo _config;

    BoundaryMatchAlgo* matchalgo_tight;
    BoundaryMatchAlgo* matchalgo_loose;
    
    BMTrackCluster3D process2Dtrack( std::vector< larlitecv::BMTrackCluster2D >& track2d, 
				     const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v );
    bool compare2Dtrack( const std::vector< BMTrackCluster2D >& track2d, const BMTrackCluster3D& track3d, const larcv::ImageMeta& meta,
			 float path_radius_cm, float endpt_radius_cm );

    void CollectCandidateBoundaryPixels( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchs,
      std::vector< dbscan::dbPoints >& combo_points, std::vector< std::vector< std::vector<int> > >& combo_cols,
      std::vector< larcv::Image2D >& matchedspacepts );

    void ClusterBoundaryHitsIntoEndpointCandidates( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchs, 
      const std::vector< dbscan::dbPoints >& combo_points, const std::vector< std::vector< std::vector<int> > >& combo_cols, 
      std::vector< std::vector<BoundaryEndPt> >& end_points, std::vector< larcv::Image2D>& matchedpixels );
    std::vector< BoundaryEndPt > DefineEndpointFromBoundaryCluster(  const BoundaryMuonTaggerAlgo::Crossings_t crossing_type, const dbscan::dbCluster& detspace_cluster, 
      const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchs,      
      const dbscan::dbPoints& combo_points, const std::vector< std::vector<int> >& combo_cols, std::vector<larcv::Image2D>& matchedpixels );    

 
    void GenerateEndPointMetaData( const std::vector< std::vector< BoundaryEndPt > >& endpts, const std::vector<larcv::Image2D>& img_v,
      const int rmax_window, const int rmin_window, const int col_window, 
      std::vector< std::vector<ClusterExtrema_t> >& candidate_metadata );
    void GetClusterExtrema( const dbscan::dbCluster& cluster, const dbscan::dbPoints& hits, ClusterExtrema_t& extrema );
    void SelectOnlyTrackEnds( const std::vector< std::vector< BoundaryEndPt > >& endpts, 
      const std::vector<larcv::Image2D>& img_v, const int rmax_window, const int rmin_window, const int col_width,
      std::vector< int >& cluster_passed );
    void CheckClusterExtrema(  const std::vector< std::vector< BoundaryEndPt > >& endpts, 
      const std::vector< std::vector<ClusterExtrema_t> >& candidate_metadata, const std::vector<larcv::Image2D>& img_v,
      std::vector<int>& passes_check );


  };


};

#endif
