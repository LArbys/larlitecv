#ifndef __CLUSTERGROUP_H__
#define __CLUSTERGROUP_H__

#include <vector>

// LArCV
#include "DataFormat/Image2D.h"
#include "DataFormat/Pixel2DCluster.h"
#include "DataFormat/Pixel2D.h"
#include "Base/PSet.h"

// LArCV/app
#include "dbscan/DBSCANAlgo.h"

// OpenCV
#ifndef __CINT__
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif
#endif

namespace larlitecv {

  // Data Product Representing a group of clusters

  class ClusterLink {
  public:
    ClusterLink( const larcv::Pixel2DCluster& a, const larcv::Pixel2DCluster& b, 
    	const larcv::Pixel2D& pix_a, const larcv::Pixel2D& pix_b, const float link_dist  ) 
      : first(&a), second(&b), firstpixel(pix_a), secondpixel(pix_b), dist(link_dist)
      {};
    virtual ~ClusterLink() {};
    const larcv::Pixel2DCluster* first;
    const larcv::Pixel2DCluster* second;
    const larcv::Pixel2D firstpixel;
    const larcv::Pixel2D secondpixel;
    const float dist;
  };

  class ClusterGroup {
  public:
    // default constructor
    ClusterGroup() {};
    virtual ~ClusterGroup() {};

    std::vector<int> m_cluster_indices;
    std::vector<larcv::Pixel2DCluster> m_clusters_v;   //< list of clusters in the group
    std::vector< dbscan::ClusterExtrema > m_extrema_v; //< list of cluster extrema    
    std::vector< ClusterLink > m_links_v;              //< list of links between the cluster groups (links extrema)
    float tick_width; //< tick span
    float tick_start;
    float tick_end;

  };

  typedef std::vector<ClusterGroup> PlaneClusterGroups;

  class ClusterGroupAlgoConfig {
  public:
    ClusterGroupAlgoConfig() {
      setdefaults();
    };
    virtual ~ClusterGroupAlgoConfig() {};

    std::vector<float> pixel_thresholds;
    int dbscan_cluster_minpoints;
    float dbscan_cluster_radius;
    float alldir_max_link_dist;
    float max_link_distance;
    float min_link_cosine;
    int single_cluster_group_min_npoints;

    void setdefaults();

    static ClusterGroupAlgoConfig MakeClusterGroupAlgoConfigFromPSet( const larcv::PSet& ps );

  };

  class ClusterGroupAlgo {
  public:
    ClusterGroupAlgo( const ClusterGroupAlgoConfig& config ) 
     : m_config(config) {};
    virtual ~ClusterGroupAlgo() {};

    std::vector<std::vector<ClusterGroup> > MakeClusterGroups( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const std::vector<larcv::Image2D>& tagged_v );


  protected:

    const ClusterGroupAlgoConfig m_config;

    // hold data while we work through the algo
    struct AlgoData_t {
      std::vector<larcv::Image2D>        untagged_v;
      std::vector<dbscan::dbPoints>      pixels_v;
      std::vector<dbscan::dbscanOutput>  dbscan_output_v;
      std::vector< std::vector<larcv::Pixel2DCluster> >  plane_clusters_v;
      std::vector< std::vector<dbscan::ClusterExtrema> > plane_extrema_v;
      std::vector< std::vector<ClusterLink> >            plane_links_v;
      std::vector< std::vector<ClusterGroup > >          plane_groups_v;
      std::vector< ClusterLink* >                  GetClusterLinks( int plane, const larcv::Pixel2DCluster* pixcluster );
    };
    AlgoData_t m_stored_algodata;

    void MakeUntaggedImages( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const std::vector<larcv::Image2D>& tagged_v, AlgoData_t& data );

    void ExtractClustersFromImage( AlgoData_t& data );

    void FindClusterLinks( const std::vector<larcv::Image2D>& img_v, AlgoData_t& data );

    void AssembleClusterGroups( AlgoData_t& data );

    void getNextLinkedCluster( AlgoData_t& data, const int& plane, std::vector< const larcv::Pixel2DCluster* >& cluster_history, 
      std::set< const larcv::Pixel2DCluster* >& clustergroup, std::vector< ClusterLink* >& used_links );     

#ifndef __CINT__
#ifdef USE_OPENCV
    std::vector<cv::Mat> makeBaseClusterImageOCV( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const std::vector<larcv::Image2D>& tagged_v, const AlgoData_t* data=nullptr );
#endif
#endif

  };

}

#endif