#include "ClusterGroupAlgo.h"

#include "TRandom.h"

namespace larlitecv {

  void ClusterGroupAlgoConfig::setdefaults() {
    pixel_thresholds.resize(3,10.0);
    dbscan_cluster_minpoints = 20;
    dbscan_cluster_radius    = 10.0;
    alldir_max_link_dist     = 30.0;
    max_link_distance        = 300.0;
    min_link_cosine          = 0.9;
    single_cluster_group_min_npoints = 100;
    save_jpg_images          = false;
  }

  ClusterGroupAlgoConfig ClusterGroupAlgoConfig::MakeClusterGroupAlgoConfigFromPSet( const larcv::PSet& ps ) {
    ClusterGroupAlgoConfig cfg;
    cfg.pixel_thresholds         = ps.get<std::vector<float> >("PixelThresholds");
    cfg.dbscan_cluster_minpoints = ps.get<int>("DBScanClusterMinPoints");
    cfg.dbscan_cluster_radius    = ps.get<int>("DBScanClusterRadius");
    cfg.alldir_max_link_dist     = ps.get<float>("AllDirMaxLinkDist");
    cfg.max_link_distance        = ps.get<float>("MaxLinkDistance");
    cfg.min_link_cosine          = ps.get<float>("MinLinkCosine");
    cfg.single_cluster_group_min_npoints = ps.get<int>("SingleClusterGroupMinNumPoints");
    cfg.save_jpg_images          = ps.get<bool>("SaveJPGImages");
  }

  std::vector< std::vector<ClusterGroup> > ClusterGroupAlgo::MakeClusterGroups( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
    const std::vector<larcv::Image2D>& tagged_v ) {

    AlgoData_t data;

    MakeUntaggedImages( img_v, badch_v, tagged_v, data );
    if ( data.untagged_v.size()!=img_v.size() )
      throw std::runtime_error("ClusterGroupAlgo::MakeClusterGroups[error]. Untagged images made not equal to number of input images.");

    ExtractClustersFromImage( data );
    std::cout << "ClusterGroupAlgo::MakeClusterGroups: number of clusters" << std::endl;
    for (size_t p=0; p<img_v.size(); p++) {
      std::cout << " plane " << p << ": " << data.plane_clusters_v.at(p).size() << std::endl;
      // for (size_t idx=0; idx<data.plane_clusters_v.at(p).size(); idx++) {
      //   if ( data.plane_clusters_v.at(p).at(idx).size()<m_config.dbscan_cluster_minpoints )
      //     continue;
      //   std::cout << "    " << idx << ": " << &data.plane_clusters_v.at(p).at(idx) << std::endl;
      // }
    }

    FindClusterLinks( img_v, data );
    std::cout << "ClusterGroupAlgo::MakeClusterGroups: number of links" << std::endl;
    for (size_t p=0; p<img_v.size(); p++) {
      std::cout << " plane " << p << ": " << data.plane_links_v.at(p).size() << std::endl;
      // for (size_t idx=0; idx<data.plane_links_v.at(p).size(); idx++) {
      //   const ClusterLink& link = data.plane_links_v.at(p).at(idx);
      //   std::cout << "    " << idx << " " << &link << " [" << link.first << " <--> " << link.second << "]" << std::endl;
      // }
    }    

    AssembleClusterGroups( data );
    std::cout << "ClusterGroupAlgo::MakeClusterGroups: number of groups" << std::endl;
    for (size_t p=0; p<img_v.size(); p++) {
      const larcv::ImageMeta& meta = img_v.at(p).meta();
      for (size_t ig=0; ig<data.plane_groups_v.at(p).size(); ig++) {
        const ClusterGroup& group = data.plane_groups_v.at(p).at(ig);
        std::cout << "Plane " << p << " Group #" << ig << ", " << group.m_cluster_indices.size() << " clusters: ";
        for ( auto& idx : group.m_cluster_indices ) {
           std::cout << idx << " ";
        }
        std::cout << " range: [" << meta.pos_y(group.tick_start) << "," << meta.pos_y(group.tick_end) << "] width=" << group.tick_width;
        std::cout << std::endl;
      }
    }

    std::vector< std::vector<ClusterGroup> > groups;

    if ( m_config.save_jpg_images ) {
      std::vector<cv::Mat> cvimgs_v = makeBaseClusterImageOCV( img_v, badch_v, tagged_v, &data );
      for ( size_t p=0; p<cvimgs_v.size(); p++) {
        cv::Mat& cvimg = cvimgs_v.at(p);
        std::stringstream ss;
        if ( m_jpg_stem=="" )
          ss << "clustergroup_p" << p << ".jpg";
        else
          ss << m_jpg_stem + "_p" << p << ".jpg";
        cv::imwrite( ss.str(), cvimg );
      }
    }

    //std::swap( groups, data.plane_groups_v );
    std::swap( m_stored_algodata, data );

    return m_stored_algodata.plane_groups_v;
  }

  void ClusterGroupAlgo::MakeUntaggedImages( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
    const std::vector<larcv::Image2D>& tagged_v, ClusterGroupAlgo::AlgoData_t& data ) {
    
    data.untagged_v.clear();
    for ( size_t p=0; p<img_v.size(); p++ ) {
      const larcv::ImageMeta& meta = img_v.at(p).meta();
      const larcv::Image2D& img    = img_v.at(p);
      const larcv::Image2D& tagged = tagged_v.at(p);

      larcv::Image2D untagged( meta );
      untagged.paint(0.0);

      for (size_t r=0; r<meta.rows(); r++) {
        for (size_t c=0; c<meta.cols(); c++) {
          if ( img.pixel(r,c) > m_config.pixel_thresholds[p] && tagged.pixel(r,c)<0.5 ) 
            untagged.set_pixel(r,c,255);
        }
      }

      data.untagged_v.emplace_back( std::move(untagged) );
    }

  }

  void ClusterGroupAlgo::ExtractClustersFromImage( ClusterGroupAlgo::AlgoData_t& data ) {

    dbscan::DBSCANAlgo algo;

    for (size_t p=0; p<data.untagged_v.size(); p++) {
      dbscan::dbPoints pixels       = dbscan::extractPointsFromImage( data.untagged_v.at(p), m_config.pixel_thresholds[p] );
      dbscan::dbscanOutput dboutput = algo.scan( pixels, m_config.dbscan_cluster_minpoints, m_config.dbscan_cluster_radius, false );
      std::vector<dbscan::ClusterExtrema> ex_v;
      std::vector<larcv::Pixel2DCluster>  cluster_v;
      for (size_t ic=0; ic<dboutput.clusters.size(); ic++) {

        // make a cluster object
        larcv::Pixel2DCluster cluster;
        for ( size_t ihit=0; ihit<dboutput.clusters.at(ic).size(); ihit++ ) {
          int hitidx = dboutput.clusters.at(ic).at(ihit);
          larcv::Pixel2D pix( pixels.at(hitidx)[0], pixels.at(hitidx)[1] );
          pix.Intensity( data.untagged_v.at(p).pixel( pix.Y(), pix.X() ) );
          cluster.emplace_back( std::move(pix) );
        }
        cluster_v.emplace_back( std::move(cluster) );

        dbscan::ClusterExtrema ex = dbscan::ClusterExtrema::FindClusterExtrema( ic, dboutput, pixels );
        ex_v.emplace_back( std::move(ex) );
      }
      // stuff it into the AlgoData struct
      data.pixels_v.emplace_back( std::move(pixels) );
      data.dbscan_output_v.emplace_back( std::move(dboutput) );
      data.plane_extrema_v.emplace_back( std::move(ex_v) );
      data.plane_clusters_v.emplace_back( std::move(cluster_v) );
    }//end of plane loop

  }

  void ClusterGroupAlgo::FindClusterLinks( const std::vector<larcv::Image2D>& img_v, ClusterGroupAlgo::AlgoData_t& data ) {
    // Build links between clusters
    const size_t nplanes = img_v.size();
    for (size_t p=0; p<nplanes; p++) {
        std::vector<ClusterLink> links_v;
      const dbscan::dbClusters& clusters = data.dbscan_output_v.at(p).clusters;
      for (size_t ic_a=0; ic_a<clusters.size(); ic_a++ ) {

        if ( (int)clusters.at(ic_a).size()<m_config.dbscan_cluster_minpoints ) continue;

        for (size_t ic_b=ic_a+1; ic_b<clusters.size(); ic_b++ ) {

          if ( (int)clusters.at(ic_b).size()<m_config.dbscan_cluster_minpoints ) continue;

          const dbscan::ClusterExtrema& ex_a = data.plane_extrema_v.at(p).at(ic_a);
          const dbscan::ClusterExtrema& ex_b = data.plane_extrema_v.at(p).at(ic_b);

          // we find closest distance between cluster extrema
          float min_dist = -1;
          dbscan::ClusterExtrema::Extrema_t min_exa;
          dbscan::ClusterExtrema::Extrema_t min_exb;
          float dir_min[2] = {0};
          for (int i=0; i<(int)dbscan::ClusterExtrema::kNumExtrema; i++) {
            for (int j=0; j<(int)dbscan::ClusterExtrema::kNumExtrema; j++) {

              float d = 0;
              for (int v=0; v<2; v++) {
                float dv = ex_a.extrema((dbscan::ClusterExtrema::Extrema_t)i)[v]-ex_b.extrema((dbscan::ClusterExtrema::Extrema_t)j)[v];
                d += dv*dv;
              }
              d = sqrt(d);
              if ( min_dist<0 || d<min_dist ) {
                min_dist = d;
                min_exa = (dbscan::ClusterExtrema::Extrema_t)i;
                min_exb = (dbscan::ClusterExtrema::Extrema_t)j;
                for (int v=0; v<2; v++)
                  dir_min[v] = (ex_a.extrema((dbscan::ClusterExtrema::Extrema_t)i)[v]-ex_b.extrema((dbscan::ClusterExtrema::Extrema_t)j)[v])/min_dist;
              }
            }
          }

          if ( min_dist <0 || min_dist > m_config.max_link_distance ) continue;

          // min-dist criteria passed. now check direction compatibility.
          // the direction of each cluster is made by using the largest distance from link-extrema
          float max_dist_exa = 0.;
          float max_dist_exb = 0.;
          int max_exa = 0;
          int max_exb = 0;
          float dir_a[2] = {0};
          float dir_b[2] = {0};

          for (int i=0; i<(int)dbscan::ClusterExtrema::kNumExtrema; i++) {
            float d_a = 0.;
            float d_b = 0.;
            for (int v=0; v<3; v++) {
              float dv_a = ex_a.extrema((dbscan::ClusterExtrema::Extrema_t)i)[v] - ex_a.extrema(min_exa)[v];
              d_a += dv_a*dv_a;
              float dv_b = ex_b.extrema((dbscan::ClusterExtrema::Extrema_t)i)[v] - ex_b.extrema(min_exb)[v];
              d_b += dv_b*dv_b;
            }
            d_a = sqrt(d_a);
            d_b = sqrt(d_b);
            if ( d_a>max_dist_exa ) {
              max_exa = i;
              max_dist_exa = d_a;
              for (int v=0; v<2; v++)
                dir_a[v] = (ex_a.extrema((dbscan::ClusterExtrema::Extrema_t)max_exa)[v] - ex_a.extrema(min_exa)[v])/max_dist_exa;
            }
            if ( d_b>max_dist_exb ) {
              max_exb = i;
              max_dist_exb = d_b;
              for (int v=0; v<2; v++)
                dir_b[v] = (ex_b.extrema((dbscan::ClusterExtrema::Extrema_t)max_exb)[v] - ex_b.extrema(min_exb)[v])/max_dist_exb;              
            }
          }

          float coscluster_a = 0.;
          float coscluster_b = 0.;          
          for (int v=0; v<2; v++) {
            coscluster_a +=  dir_min[v]*dir_a[v];
            coscluster_b += -dir_min[v]*dir_b[v];
          }

          // cosine must be above some value for both (or within some distance)
          if ( min_dist > m_config.alldir_max_link_dist && (coscluster_a<m_config.min_link_cosine || coscluster_b<m_config.min_link_cosine) )
            continue;

          // define the link
          larcv::Pixel2D pix_a( ex_a.extrema(min_exa)[0], ex_a.extrema(min_exa)[1] );
          larcv::Pixel2D pix_b( ex_b.extrema(min_exb)[0], ex_b.extrema(min_exb)[1] );          
          ClusterLink link( data.plane_clusters_v.at(p).at(ic_a), data.plane_clusters_v.at(p).at(ic_b), pix_a, pix_b, min_dist );

          links_v.emplace_back( std::move(link) );
        }// end of cluster_b loop
      } // end of cluster_a loop
      data.plane_links_v.emplace_back( std::move(links_v) );
    } //end loop over planes

  }

  void ClusterGroupAlgo::AssembleClusterGroups( ClusterGroupAlgo::AlgoData_t& data ) {

    for ( size_t p=0; p<data.plane_clusters_v.size(); p++ ) {
      const larcv::ImageMeta& meta = data.untagged_v.at(p).meta();
      std::vector< ClusterGroup > plane_group_v;
      std::set<int> used_cluster_indices;
      while ( true ) {

        // find a new cluster index
        int seed_index = -1;
        for ( int ic=0; ic<(int)data.plane_clusters_v.at(p).size(); ic++ ) {
          if ( data.plane_clusters_v.at(p).at(ic).size()<m_config.dbscan_cluster_minpoints )
            continue;
          if ( used_cluster_indices.find(ic)==used_cluster_indices.end() ) {
            seed_index = ic;
            break;
          }
        }
  
        if ( seed_index<0 )
          break;
  
        used_cluster_indices.insert(seed_index);
        larcv::Pixel2DCluster& pseed_cluster = data.plane_clusters_v.at(p).at(seed_index);
  
        std::vector< larcv::Pixel2DCluster const* > cluster_history;
        cluster_history.push_back( &pseed_cluster );
        std::set< larcv::Pixel2DCluster const* > cluster_group;
        cluster_group.insert( &pseed_cluster );
        std::vector< ClusterLink* > used_links;
        getNextLinkedCluster( data, p, cluster_history, cluster_group, used_links );

        if ( cluster_group.size()<=1 && (int)pseed_cluster.size()<m_config.single_cluster_group_min_npoints )
          continue;

        // make the data product, by copy (barf!)
        ClusterGroup cgroup;

        // match pointers to get the cluster and cluster extrema. this is so gross.
        for ( int idx=0; idx<(int)data.plane_clusters_v.at(p).size(); idx++) {
          larcv::Pixel2DCluster* pcluster = &(data.plane_clusters_v.at(p).at(idx));
          if ( cluster_group.find( pcluster )!=cluster_group.end() ) {
            cgroup.m_extrema_v.push_back( data.plane_extrema_v.at(p).at(idx) );
            cgroup.m_clusters_v.push_back( *pcluster );
            cgroup.m_cluster_indices.push_back(idx);
            used_cluster_indices.insert(idx);
          }
        }

        // determine time span
        cgroup.tick_start = -1;
        cgroup.tick_end   = -1;
        for ( auto const& ex : cgroup.m_extrema_v ) {
          //std::cout << "  top=" << meta.pos_y(ex.topmost()[1]) << " bot=" << meta.pos_y(ex.bottommost()[1]) << std::endl;
          if ( cgroup.tick_start<0 || cgroup.tick_start>ex.bottommost()[1] )
            cgroup.tick_start = ex.bottommost()[1];
          if ( cgroup.tick_end<0 || cgroup.tick_end<ex.topmost()[1] )
            cgroup.tick_end = ex.topmost()[1];
        }
        cgroup.tick_width = cgroup.tick_end - cgroup.tick_start;

        // pass on the links
        for ( auto & plink : used_links ) {
          cgroup.m_links_v.push_back( *plink );
        }

        // store the group
        plane_group_v.emplace_back( std::move(cgroup) );
      } 

      // store groups found on plane
      data.plane_groups_v.emplace_back( std::move(plane_group_v) );
    }//end of plane loop

  }

  void ClusterGroupAlgo::getNextLinkedCluster( ClusterGroupAlgo::AlgoData_t& data, const int& plane, 
    std::vector< const larcv::Pixel2DCluster* >& cluster_history, std::set< const larcv::Pixel2DCluster* >& clustergroup, 
    std::vector< ClusterLink* >& used_links ) {

    // get the links
    if ( cluster_history.size()==0)
      return;
    const larcv::Pixel2DCluster* current_cluster = cluster_history.back();
    std::vector< ClusterLink*> cluster_link = data.GetClusterLinks( plane, current_cluster );
    for ( auto &plink : cluster_link ) {
      // go to first cluster not already in the group
      if ( plink->first==current_cluster && clustergroup.find( plink->second )==clustergroup.end() ) {
        clustergroup.insert( plink->second );
        cluster_history.push_back( plink->second );
        used_links.push_back( plink );
        getNextLinkedCluster( data, plane, cluster_history, clustergroup, used_links );
      }
      else if ( plink->second==current_cluster && clustergroup.find( plink->first )==clustergroup.end() ) {
        clustergroup.insert( plink->first );
        cluster_history.push_back( plink->first );
        used_links.push_back( plink );
        getNextLinkedCluster( data, plane, cluster_history, clustergroup, used_links );
      }
    }
    // no link
    cluster_history.pop_back();
    return;
  }  

  std::vector< ClusterLink* > ClusterGroupAlgo::AlgoData_t::GetClusterLinks( int plane, const larcv::Pixel2DCluster* pixcluster ) {
    std::vector< ClusterLink* > links_v;
    // we simply cycle through the list of clusters and get those links attached to this cluster
    for ( auto & link : plane_links_v.at(plane) ) {
      if ( pixcluster==link.first || pixcluster==link.second ) {
        links_v.push_back( &link );
      }
    }
    return links_v;
  }

#ifndef __CINT__
#ifdef USE_OPENCV
  std::vector<cv::Mat> ClusterGroupAlgo::makeBaseClusterImageOCV( const std::vector<larcv::Image2D>& img_v, 
    const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v, const AlgoData_t* data  ) {

    if ( data==nullptr ) {
      data = &m_stored_algodata;
    }

    // we draw an image, that highlights the clusters, links between them, and interesting space points
    // we drawn an image per plane

    TRandom rand(1);
    std::vector<cv::Mat> cvimgs_v;

    for ( size_t p=0; p<img_v.size(); p++ ) {

      const larcv::Image2D& img = img_v.at(p);
      const larcv::ImageMeta& meta = img.meta();

      // first make a CV image we can have fun with
      cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img, m_config.pixel_thresholds[p], 100 );
      
      // color in bad pixels
      const larcv::Image2D& bad = badch_v.at(p);
      for (size_t r=0; r<bad.meta().rows(); r++) {
        for (size_t c=0; c<bad.meta().cols(); c++) {
          if ( bad.pixel(r,c)>0 ){
            cv::Vec3b& pix = cvimg.at<cv::Vec3b>( cv::Point(c,r) );
            pix[0] = 50;
            pix[1] = 50;
            pix[2] = 50;
          }
        }
      }

      // color in tagged pixels
      const larcv::Image2D& tagged = tagged_v.at(p);
      for (size_t r=0; r<tagged.meta().rows(); r++) {
        for (size_t c=0; c<tagged.meta().cols(); c++) {
          if ( tagged.pixel(r,c)>0 ){
            cv::Vec3b& pix = cvimg.at<cv::Vec3b>( cv::Point(c,r) );
            pix[0] = 255;
            pix[1] = 0;
            pix[2] = 0;
          }
        }
      }

      // draw the links
      for ( auto const& link : data->plane_links_v.at(p) ) {
        cv::line( cvimg, cv::Point(link.firstpixel.X(),link.firstpixel.Y()), cv::Point(link.secondpixel.X(),link.secondpixel.Y()), cv::Scalar(0,0,150), 2 );
      }

      // ok, we now label base clusters with colors! also draw their extrema and their links     
      const std::vector< larcv::Pixel2DCluster >& clusters = data->plane_clusters_v.at(p);
      const std::vector< dbscan::ClusterExtrema >& extrema = data->plane_extrema_v.at(p);

      for (size_t icluster=0; icluster<clusters.size(); icluster++){
        if ( clusters.at(icluster).size()<m_config.dbscan_cluster_minpoints ) continue;
        const larcv::Pixel2DCluster& cluster = clusters.at(icluster);
        // pick a color
        cv::Vec3b color;
        color[0] = (int)(rand.Uniform()*255);
        color[1] = (int)(rand.Uniform()*255);
        color[2] = (int)(rand.Uniform()*255);
        for (auto const& pix : cluster ) {
          cvimg.at<cv::Vec3b>( cv::Point(pix.X(),pix.Y()) ) = color;
        }
        // color in the extrema
        const dbscan::ClusterExtrema& ex = extrema.at(icluster);
        cv::circle(cvimg, cv::Point(ex.leftmost()[0],   ex.leftmost()[1]),   3, cv::Scalar(0,255,0),-1);
        cv::circle(cvimg, cv::Point(ex.topmost()[0],    ex.topmost()[1]),    3, cv::Scalar(255,255,0),-1);          
        cv::circle(cvimg, cv::Point(ex.rightmost()[0],  ex.rightmost()[1]),  3, cv::Scalar(0,255,255),-1);          
        cv::circle(cvimg, cv::Point(ex.bottommost()[0], ex.bottommost()[1]), 3, cv::Scalar(255,0,255),-1);

        std::stringstream slabel;
        slabel << icluster;
        cv::putText(cvimg,slabel.str(),cv::Point(ex.leftmost()[0],ex.leftmost()[1]),cv::FONT_HERSHEY_PLAIN, 1, cv::Scalar(255,255,255) );
      }//end of cluster loop

      cvimgs_v.emplace_back( std::move(cvimg) );
    }

    return cvimgs_v;
  }  
#endif
#endif

}
