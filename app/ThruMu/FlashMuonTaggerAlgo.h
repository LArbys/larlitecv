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

// larlitecv header
#include "BoundaryMuonTaggerTypes.h"
#include "BoundaryEndPt.h"
#include "BoundarySpacePoint.h"
#include "FlashMuonTaggerConfig.h"

namespace larlitecv {
    
  class FlashMuonTaggerAlgo {

  public:
    
    typedef enum { kAnode, kCathode, kOutOfImage } SearchMode_t;
    typedef enum { kFront, kBack, kNotSet } OutOfImage_t;

    FlashMuonTaggerAlgo( SearchMode_t mode_) {
      fSearchMode = mode_;
      fOutOfImageMode = kNotSet;
      fconfigured = false;
      loadGeoInfo();
    };
    virtual ~FlashMuonTaggerAlgo() {};
    
    struct ClusterInfo_t {
      // this needs to inherit from ClusterExtrema at some point.
      int plane; //< plane it was found on
      int cluster_idx; //< index of the cluster in the dbscan output
      int hits_rmin_boundary;
      int hits_rmax_boundary;
      int ntime_ends; //< how many boundnaries it reached
      int row;
      int col;
      float q;
      bool indead_region;
      int row_min;
      int col_min;
      float q_min;
      int row_max;
      int col_max;
      float q_max;
      int npixels;
    };    

    void configure( const FlashMuonTaggerConfig& config_ ) { fConfig=config_; fconfigured = true; };

    bool flashMatchTrackEnds( const std::vector< larlite::event_opflash* >& opflashsets, const std::vector<larcv::Image2D>& tpc_imgs,
			      const std::vector<larcv::Image2D>& badch_imgs,
			      std::vector< BoundarySpacePoint >& trackendpts );

    bool findImageTrackEnds( const std::vector<larcv::Image2D>& tpc_imgs, const std::vector<larcv::Image2D>& badch_imgs,
			     std::vector< BoundarySpacePoint >& trackendpts );

    BoundaryEnd_t SearchModeToEndType( SearchMode_t mode );
    
  protected:
    
    SearchMode_t fSearchMode;
    OutOfImage_t fOutOfImageMode;
    bool fconfigured;
    FlashMuonTaggerConfig fConfig;
    
    // subroutines
    bool findClusterEnds( const dbscan::dbscanOutput& clout, const dbscan::dbPoints& winpoints, 
			  const int clusterid, const int row_target, const int plane, const larcv::ImageMeta& meta,
			  BoundaryEndPt& endpt );

    std::vector<FlashMuonTaggerAlgo::ClusterInfo_t> analyzeClusters( const dbscan::dbscanOutput& dbscan_output, 
      const dbscan::dbPoints& hits, const larcv::Image2D& img, const int row_target, const int plane, const int row_radius );

    void findPlaneMatchedClusters( const std::vector< std::vector< FlashMuonTaggerAlgo::ClusterInfo_t > >& cluster_info, 
      const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
      const float max_triarea, const std::vector<float>& z_range, const std::vector<float>& y_range, std::vector< BoundarySpacePoint >& endpts_v,
      std::vector< std::vector< FlashMuonTaggerAlgo::ClusterInfo_t > >& accepted_clusters );

    void generate3PlaneIntersections(  const std::vector< std::vector< ClusterInfo_t > >& cluster_info, const std::vector<larcv::Image2D>& img_v, 
      const std::vector<float>& z_range, const std::vector<float>& y_range, const float max_triarea,
      std::vector< BoundarySpacePoint >& endpts_v, std::vector< std::vector< ClusterInfo_t > >& accepted_clusters,  
      std::vector< std::vector<int> >& cluster_used );

    void generate2PlaneIntersections( const std::vector< std::vector< ClusterInfo_t > >& cluster_info, 
      const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
      const std::vector<float>& z_range, const std::vector<float>& y_range,
      std::vector< BoundarySpacePoint >& endpts_v, std::vector< std::vector< ClusterInfo_t > >& accepted_clusters,
      std::vector< std::vector<int> >& cluster_used );


    void filterClusters( const std::vector< std::vector< ClusterInfo_t > >& accepted_cluster_matches, 
      const std::vector<larcv::Image2D>& img_v, const int rmax_window, const int rmin_window, const int col_width, std::vector< int >& cluster_passed );    


    void loadGeoInfo();

    
    void FindFlashesByChargeClusters( const int row_target, const larlitecv::BoundaryEnd_t point_type,
				      const std::vector<larcv::Image2D>& tpc_imgs, const std::vector<larcv::Image2D>& badch_imgs,
				      const std::vector<float>& z_range, const std::vector<float>& y_range, std::vector< BoundarySpacePoint >& trackendpts );

    void FindFlashesBy3DSegments( const int row_target, const larlitecv::BoundaryEnd_t point_type,
				  const std::vector<larcv::Image2D>& tpc_imgs, const std::vector<larcv::Image2D>& badch_imgs,
				  const std::vector<float>& z_range, std::vector< BoundarySpacePoint >& trackendpts );
    
    

  public:
    int fNPMTs;
    float pmtpos[32][3];    
    
  };
  
  
}

#endif
