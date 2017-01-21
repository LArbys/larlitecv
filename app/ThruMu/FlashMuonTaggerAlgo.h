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
    int                 search_row_radius; ///< neighborhood to gather charge for clustering
    float               flash_zrange_extension;
    int                 flash_pixelcluster_minsize; //< min. number of pixels required to make a flash-tagged cluster
    float               max_triarea;
    float               max_triarea_tight;
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
    
    struct ClusterInfo_t {
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
			      std::vector< std::vector< BoundaryEndPt > >& trackendpts );
    bool findImageTrackEnds( const std::vector<larcv::Image2D>& tpc_imgs, const std::vector<larcv::Image2D>& badch_imgs,
			     std::vector< std::vector< BoundaryEndPt > >& trackendpts );
    
  protected:
    
    SearchMode_t fSearchMode;
    bool fconfigured;
    FlashMuonTaggerConfig fConfig;
    
    // subroutines
    bool findClusterEnds( const dbscan::dbscanOutput& clout, const dbscan::dbPoints& winpoints, 
			  const int clusterid, const int row_target, const int plane, const larcv::ImageMeta& meta,
			  BoundaryEndPt& endpt );
    void defineClusterPositions( const dbscan::dbscanOutput& clout, const dbscan::dbPoints& winpoints, const int row_target, 
                                 const int plane, const int point_type, const larcv::Image2D& img, std::vector<BoundaryEndPt>& endpt_v );
    std::vector< std::vector<BoundaryEndPt> > generate3PlaneIntersections( const std::vector< std::vector<BoundaryEndPt> >& endptset_v, 
      const float max_triarea, const std::vector<larcv::Image2D>& img_v, const std::vector<float>& z_range, const std::vector<float>& y_range );
    std::vector<FlashMuonTaggerAlgo::ClusterInfo_t> analyzeClusters( const dbscan::dbscanOutput& dbscan_output, 
      const dbscan::dbPoints& hits, const larcv::Image2D& img, const int row_target, const int plane, const int row_radius );
    void findPlaneMatchedClusters( const std::vector< std::vector< FlashMuonTaggerAlgo::ClusterInfo_t > >& cluster_info, const std::vector<larcv::Image2D>& img_v,
      const float max_triarea, const std::vector<float>& z_range, const std::vector<float>& y_range, std::vector< std::vector<BoundaryEndPt> >& endpts_v,
      std::vector< std::vector< FlashMuonTaggerAlgo::ClusterInfo_t > >& accepted_clusters );

    void filterClusters( const std::vector< std::vector< ClusterInfo_t > >& accepted_cluster_matches, 
      const std::vector<larcv::Image2D>& img_v, const int rmax_window, const int rmin_window, const int col_width, std::vector< int >& cluster_passed );


    void loadGeoInfo();

    
  public:
    int fNPMTs;
    float pmtpos[32][3];    
    
  };
  
  
}

#endif
