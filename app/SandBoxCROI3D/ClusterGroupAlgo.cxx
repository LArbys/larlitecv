#include "ClusterGroup.h"

namespace larlitecv {

  void ClusterGroupAlgoConfig::setdefaults() {
    pixel_thresholds.resize(3,10.0);
    dbscan_cluster_minpoints = 10;
    dbscan_cluster_radius    = 10.0;
    alldir_max_link_dist     = 20.0;
    max_link_distance        = 100.0;
    min_link_cosine          = 0.75;
  }


  std::vector<ClusterGroup> ClusterGroupAlgo::MakeClusterGroups( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
    const std::vector<larcv::Image2D>& tagged_v ) {

    AlgoData_t data;

    MakeUntaggedImages( img_v, badch_v, tagged_v, data );

    ExtractClustersFromImage( data );

    FindClusterLinks( img_v, data );

    AssembleClusterGroups( data );

    std::vector<ClusterGroup> groups;

    return groups;
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

          if ( min_dist>m_config.alldir_max_link_dist && (min_dist <0 || min_dist > m_config.max_link_distance ) ) continue;

          // min-dist criteria passed. now check direction compatibility
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

          // cosine must be above some value for both
          if ( coscluster_b<m_config.min_link_cosine )
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

      std::set<int> used_cluster_indices;
      while ( true ) {

        // find a new cluster index
        int seed_index = -1;
        for ( int ic=0; ic<(int)data.plane_clusters_v.size(); ic++ ) {
          if ( used_cluster_indices.find(ic)==used_cluster_indices.end() ) {
            seed_index = ic;
            break;
          }
        }
  
        if ( seed_index<0 )
          break;
  
        used_cluster_indices.insert(seed_index);
        const larcv::Pixel2DCluster* pseed_cluster = &( data.plane_clusters_v.at(p).at(seed_index) );
  
        ClusterGroup cgroup;
        std::vector< const larcv::Pixel2DCluster* > cluster_history;
        cluster_history.push_back( pseed_cluster );
        std::set< const larcv::Pixel2DCluster* > cluster_group;
        std::vector< const ClusterLink* > used_links;
        getNextLinkedCluster( data, p, cluster_history, cluster_group, used_links );
      }	
    }//end of plane loop

  }

  void ClusterGroupAlgo::getNextLinkedCluster( ClusterGroupAlgo::AlgoData_t& data, const int& plane, 
    std::vector<const larcv::Pixel2DCluster* >& cluster_history, std::set<const larcv::Pixel2DCluster* >& clustergroup, 
    std::vector< const ClusterLink* >& used_links ) {

    // get the links
    if ( cluster_history.size()==0)
      return;
    const larcv::Pixel2DCluster* current_cluster = cluster_history.back();
    std::vector<const ClusterLink*> cluster_link = data.GetClusterLinks( plane, current_cluster );
    for ( auto const& plink : cluster_link ) {
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

  std::vector< const ClusterLink* > ClusterGroupAlgo::AlgoData_t::GetClusterLinks( int plane, const larcv::Pixel2DCluster* pixcluster ) {
    std::vector< const ClusterLink* > links_v;
    // we simply cycle through the list of clusters and get those links attached to this cluster
    for ( auto const& link : plane_links_v.at(plane) ) {
      if ( pixcluster==link.first || pixcluster==link.second ) {
        links_v.push_back( &link );
      }
    }
    return links_v;
  }

}