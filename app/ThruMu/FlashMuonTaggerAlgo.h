#ifndef __FLASH_MUON_TAGGER_ALGO__
#define __FLASH_MUON_TAGGER_ALGO__

#include <vector>

// larlite data
#include "Base/DataFormatConstants.h"
#include "DataFormat/opflash.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "ANN/ANNAlgo.h"
#include "dbscan/DBSCANAlgo.h"
#include "PMTWeights/WireData.h"

// larlitecv header
#include "BoundaryEndPt.h"

namespace larlitecv {
  
  class FlashMuonTaggerConfig {
  public:
    FlashMuonTaggerConfig() {};
    virtual ~FlashMuonTaggerConfig() {};
    
    // values per plane
    std::vector<float>  pixel_value_threshold;
    std::vector<int>    clustering_time_neighborhood;
    std::vector<int>    clustering_wire_neighborhood;
    std::vector<int>    clustering_minpoints;
    std::vector<double> clustering_radius;
    std::vector<int>    endpoint_time_neighborhood;
    int                 verbosity;
    float               trigger_tick;
    float               usec_per_tick;
    float               drift_distance;
    float               drift_velocity;
  };
  
  class FlashMuonTaggerAlgo {

  public:
    
    typedef enum { kAnode, kCathode, kOutOfImage } SearchMode_t;

    FlashMuonTaggerAlgo( SearchMode_t mode_) {
      fSearchMode = mode_;
      fconfigured = false;
      loadGeoInfo();
    };
    virtual ~FlashMuonTaggerAlgo() {};
    
    void configure( const FlashMuonTaggerConfig& config_ ) { fConfig=config_; fconfigured = true; };
    bool findTrackEnds( const std::vector< larlite::event_opflash* >& opflashsets, const larcv::Image2D& image,
			std::vector< BoundaryEndPt >& trackendpts, larcv::Image2D& markedimg );
    bool flashMatchTrackEnds( const std::vector< larlite::event_opflash* >& opflashsets, const std::vector<larcv::Image2D>& tpc_imgs,
			      std::vector< std::vector< BoundaryEndPt > >& trackendpts, std::vector< larcv::Image2D >& markedimg );
    bool findImageBoundaryEnds( const larcv::Image2D& tpc_img, std::vector< BoundaryEndPt >& trackendpts, larcv::Image2D& markedimg );    
    bool findImageTrackEnds( const std::vector<larcv::Image2D>& tpc_imgs, std::vector< std::vector< BoundaryEndPt > >& trackendpts, std::vector< larcv::Image2D >& markedimgs );
    
  protected:
    
    SearchMode_t fSearchMode;
    bool fconfigured;
    FlashMuonTaggerConfig fConfig;
    
    // subroutines
    bool findClusterEnds( const dbscan::dbscanOutput& clout, const dbscan::dbPoints& winpoints, 
			  const int clusterid, const int row_target, const int plane, const larcv::ImageMeta& meta,
			  BoundaryEndPt& endpt, larcv::Image2D& markedimg  );
    bool getClusters( const larcv::Image2D& tpc_img, int query_row, int query_col, dbscan::dbscanOutput& cluster_info, dbscan::dbPoints& windowpts, int& containing_cluster );
    void loadGeoInfo();
    float calculateIntersectionTriangle( const std::vector<float>& p0, const std::vector<float>& p1, const std::vector<float>& p2 );
    void findWireIntersections( const std::vector< std::vector<int> >& wirelists, 
				const std::vector< std::vector<float> >& valid_range,
				std::vector< std::vector<int> >& intersections3plane,
				std::vector< std::vector<float> >& vertex3plane,
				std::vector<float>& areas3plane,
				std::vector< std::vector<int> >& intersections2plane, std::vector< std::vector<float> >& vertex2plane );

  public:
    void lineSegmentIntersection2D( const std::vector< std::vector<float> >& ls1, const std::vector< std::vector<float> >& ls2, std::vector<float>& insec, int& crosses );
    
  public:
    // geometric data
    std::map<int,larcv::pmtweights::WireData> m_WireData; // key is plane ID, value is class with wire info
    int fNPMTs;
    float pmtpos[32][3];    
    
  };
  
  
}

#endif
