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
#include "ThruMu/BMTrackCluster3D.h"
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

    void findStopMuTracks( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const std::vector<larcv::Image2D>& thrumu_v, const std::vector< std::vector< const larcv::Pixel2D* > >& endpts_v );


    // We run the algorithm in passes. Being less restrictive each time.
    // Each pass has data we collect
    typedef std::vector< std::vector<AStar3DNode> > PathList_t;    
    struct PassOutput_t {
      std::vector< std::vector<larcv::Image2D> > m_cluster_images; // (deprecate)
      std::vector<larcv::Image2D> m_masked_v; // image where pixels are tagged by thrumu or previous stopmu passes
      std::vector<untagged_cluster_info_t> m_untagged_clusters_v; // clusters from dbscan over m_masked_v
      std::vector<ClusterGroup_t> m_clustergroups;  //  secondary clustering of untagged_clusters_v clusters
      std::vector<BoundarySpacePoint> m_spacepoints; // space points found by finding 3D-consistent points in the cluster groups
      PathList_t m_paths; // list of astar paths
      std::vector<int> m_path_goalreached; // marks those that reach their goal
      std::vector<int> m_endpt_index; // index of endpt used to start the path
    };
    std::vector< PassOutput_t > m_pass_data;

    // ------------------------------------
    // PASS METHODS

    PassOutput_t performPass( const StopMuClusterConfig::PassConfig_t& passcfg, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const std::vector<larcv::Image2D>& thrumu_v, const std::vector< std::vector< const larcv::Pixel2D* > >& endpts_v, std::vector<int>& endpts_used );

    // pass sub-methods

    void extractBaseClusters( const StopMuClusterConfig::PassConfig_t& passcfg, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v, 
      const std::vector< std::vector< const larcv::Pixel2D* > >& endpts, PassOutput_t& output );

    void findClusterLinks( const StopMuClusterConfig::PassConfig_t& passcfg, const std::vector<larcv::Image2D>& img_v, PassOutput_t& data );

    void analyzeClusters( const StopMuClusterConfig::PassConfig_t& passcfg, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const std::vector<larcv::Image2D>& thrumu_v, const std::vector< std::vector< const larcv::Pixel2D* > >& endpts_v, std::vector<int>& endpts_used, 
      StopMuCluster::PassOutput_t& data);

    std::vector<BoundarySpacePoint> generateCluster2PlaneSpacepoints( const StopMuClusterConfig::PassConfig_t& passcfg, const ClusterGroup_t& cluster_group, 
      const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& thrumu_v, PassOutput_t& data );

    std::vector<BoundarySpacePoint> generateCluster3PlaneSpacepoints( const StopMuClusterConfig::PassConfig_t& passcfg, const std::vector<larcv::Image2D>& img_v,
      const std::vector< std::set<int> >& cluster_groups, PassOutput_t& data );

    void saveClusterImageOCV( std::string filename ); ///< dumps out image of intermediate quantities in algorithm   

    // ------------------------------------------------------------------

    
    void runAStar();

    void postProcessing();

    void packageOutput(); 

    void getNextLinkedCluster( StopMuCluster::PassOutput_t& data, const int& plane, std::vector<int>& cluster_history, 
      std::set<int>& clustergroup, std::vector< const ClusterLink_t* >& used_links );

    BMTrackCluster3D makeBMTrackCluster3D( const std::vector<AStar3DNode>& path, 
      const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<const larcv::Pixel2D*>& start_pt );

    int m_verbosity;
    StopMuClusterConfig m_config;
    void setVerbosity( int v ) { m_verbosity = v; };    

#ifdef USE_OPENCV
    std::vector<cv::Mat> makeBaseClusterImageOCV( const PassOutput_t& data, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v );
#endif

  protected:
    /*
    float _norm( std::vector<float>& vec );
    void _wire2pixel( const int tick, const std::vector<int>& wid, const larcv::ImageMeta& meta, std::vector<int>& pixel_col, int& pixel_row );
    */
  };

}


#endif
