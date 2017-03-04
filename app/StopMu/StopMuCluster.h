#ifndef __STOPMU_CLUSTER__
#define __STOPMU_CLUSTER__

/*! \brief Clustering Algo for Stopping Muons
 *         
 *  We take advantage of previous through-going muon tagging.
 *  Given original image and image from thrumu-tagger, we cluster remaining pixels using DBSCAN
 *  We then analyze the clusters, finding their extrema.
 *  We then do an ugly O( N*(N-1)*4*4 ) operation to determine potential links between clusters.
 *  We build up a list of such links between clusters.
 *
 *  Given a starting point (a boundary crossing point from previous steps), we
 *    1) ID the cluster the starting point belongs to
 *    2) build out the largest possible cluster-group using the links
 *    3) ID extrema that form consistent 3D points between the clusters
 *    4) Building 3D space points also through finding 2-plane matches and checking for charge presence in other plane
 *    x) What I really need is an end-of-track feature finder!
 */

#include <vector>
#include <set>
#include <array>
#include <exception>
#include <algorithm>
#include <utility>

// OpenCV
#ifndef __CINT__
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif
#endif

// larcv
#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"
#include "DataFormat/Pixel2D.h"

// larcv/app
#include "dbscan/DBSCANAlgo.h"

#include "ThruMu/BoundarySpacePoint.h"
#include "ThruMu/AStar3DAlgo.h"
#include "StopMuClusterConfig.h"
#include "SMClusterTypes.h"

namespace larlitecv {


  class StopMuCluster {

    StopMuCluster() { 
      m_verbosity = 0;
    };

  public:
      
    StopMuCluster( const StopMuClusterConfig& cfg );
    virtual ~StopMuCluster() {};

    void findStopMuTrack( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const std::vector<larcv::Image2D>& thrumu_v, const std::vector< std::vector< const larcv::Pixel2D* > >& endpts_v );


    void extractBaseClusters(const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v, const std::vector< std::vector< const larcv::Pixel2D* > >& endpts );
    void saveClusterImageOCV( std::string filename ); ///< dumps out image of intermediate quantities in algorithm
    void findClusterLinks();
    void getClusterGroupForSpacepoint() {};
    std::vector<BoundarySpacePoint> generateCluster3PlaneSpacepoints( const std::vector< std::set<int> >& cluster_groups );
    std::vector<BoundarySpacePoint> generateCluster2PlaneSpacepoints( const ClusterGroup_t& cluster_group, 
      const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& thrumu_v );

    void generateCluster2PlaneSpacepoints() {};
    void runAStar();
    void postProcessing();
    void packageOutput(); 
    void getNextLinkedCluster( const int& plane, std::vector<int>& cluster_history, std::set<int>& clustergroup, std::vector< const ClusterLink_t* >& used_links );

    int m_verbosity;
    StopMuClusterConfig m_config;
    void setVerbosity( int v ) { m_verbosity = v; };    

    // stored information
    std::vector<const larcv::Image2D*> m_img_v;
    std::vector<const larcv::Image2D*> m_thrumu_v;
    std::vector< std::vector<larcv::Image2D> > m_cluster_images;
    std::vector<larcv::Image2D> m_masked_v;
    std::vector<untagged_cluster_info_t> m_untagged_clusters_v;
    std::vector<BoundarySpacePoint> m_spacepoints;
    std::vector<ClusterGroup_t> m_clustergroups;

    typedef std::vector< std::vector<AStar3DNode> > PathList_t;
    PathList_t m_paths;
    std::vector<int> m_path_goalreached;

#ifdef USE_OPENCV
    std::vector<cv::Mat> makeBaseClusterImageOCV();
#endif

  protected:
    /*
    float _norm( std::vector<float>& vec );
    void _wire2pixel( const int tick, const std::vector<int>& wid, const larcv::ImageMeta& meta, std::vector<int>& pixel_col, int& pixel_row );
    */
  };

}


#endif
