#include "StopMuCluster.h"

#include <stdio.h>
#include <time.h>
#include <sstream>

// ROOT
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

// larcv
#include "UBWireTool/UBWireTool.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larlitecv
#include "ThruMu/BoundaryEndPt.h"

namespace larlitecv {

  StopMuCluster::StopMuCluster( const StopMuClusterConfig& cfg ) : m_config(cfg) {
    // Constructor
    setVerbosity(m_config.verbosity);
    m_cvout_stem = "";
  }

  std::vector<larlitecv::BMTrackCluster3D> StopMuCluster::findStopMuTracks( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const std::vector<larcv::Image2D>& thrumu_v, const std::vector< std::vector< const larcv::Pixel2D* > >& endpts_v ) {

    std::vector<BMTrackCluster3D> stopmu_tracks;

    // make a copy of thrumu_v images -- this is the marked up images
    std::vector<larcv::Image2D> marked_v;
    for ( auto const& thrumu : thrumu_v ) {
      larcv::Image2D marked( thrumu );
      marked_v.emplace_back( std::move(marked) );
    }

    std::vector<int> endpts_used( endpts_v.size(), 0 );

    for (int ipass=0; ipass<m_config.num_passes; ipass++ ) {

      // update the marked image with stopmu_tracks
      std::cout << "===============================================================" << std::endl;
      std::cout << "[ START OF PASS " << ipass+1 << " ]" << std::endl;

      PassOutput_t pass_data = performPass( m_config.pass_configs.at(ipass), img_v, badch_v, marked_v, endpts_v, endpts_used );

      // dump out pass data
      if ( m_config.save_pass_images ) {
        std::vector<cv::Mat> passcv = makeBaseClusterImageOCV( pass_data, img_v, marked_v );
        for (size_t p=0; p<passcv.size(); p++) {
          std::stringstream ss;
          if ( m_cvout_stem=="")
            ss << "test_smc_pass" << ipass+1 << "_p" << p << ".jpg";
          else
            ss << m_cvout_stem << "_pass" << ipass+1 << "_p" << p << ".jpg";
          cv::imwrite( ss.str(), passcv.at(p) );
        }
      }

      for ( size_t ipath=0; ipath<pass_data.m_paths.size(); ipath++) {
        if ( pass_data.m_path_goalreached.at(ipath)==1 ) {
          BMTrackCluster3D track3d = makeBMTrackCluster3D( pass_data.m_paths.at(ipath), img_v, badch_v, endpts_v.at(pass_data.m_endpt_index.at(ipath)) );

          // mark up the image
          for (size_t p=0; p<img_v.size(); p++) {
            for ( auto const& pix : track3d.plane_paths.at(p).pixelpath ) {
              marked_v.at(p).set_pixel( pix.Y(), pix.X(), 255 );
            }
          }

          stopmu_tracks.emplace_back( std::move(track3d) );
        }
      }

      std::cout << "PASS " << ipass+1 << ": number of tracks=" << stopmu_tracks.size() << std::endl;
    }

    std::cout << "StopMu Tracks Found: " << stopmu_tracks.size() << std::endl;

    // dump out an image
    for (size_t p=0; p<thrumu_v.size(); p++) {
      cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img_v.at(p), 0, 50.0 );
      for (size_t c=0; c<img_v.at(p).meta().cols(); c++) {        
        for (size_t r=0; r<img_v.at(p).meta().rows(); r++) {
          bool marked = false;
          if ( marked_v.at(p).pixel(r,c)>0 ) {
            cv::Vec3b& color = cvimg.at<cv::Vec3b>( cv::Point(c,r) );
            color[0] = 0;
            color[1] = 0;
            color[2] = 255;
            marked = true;
          }
          if ( thrumu_v.at(p).pixel(r,c)>0 ) {
            cv::Vec3b& color = cvimg.at<cv::Vec3b>( cv::Point(c,r) );            
            color[0] = 255;
            color[1] = 0;
            color[2] = 0;
            marked = true;
          }
          if ( !marked && badch_v.at(p).pixel(r,c)>0 ) {
            cv::Vec3b& color = cvimg.at<cv::Vec3b>( cv::Point(c,r) );            
            color[0] = 30;
            color[1] = 30;
            color[2] = 30;
          }
        }
      }

      if ( m_config.dump_tagged_images ) {
        std::stringstream ss;
        if ( m_cvout_stem=="")
          ss << "test_smv_tagged_p" << p << ".jpg";
         else
          ss << m_cvout_stem << "_tagged_p" << p << ".jpg";
        cv::imwrite( ss.str(), cvimg );
      }

    }

    return stopmu_tracks;
  }


  StopMuCluster::PassOutput_t StopMuCluster::performPass( const StopMuClusterConfig::PassConfig_t& passcfg,  const std::vector<larcv::Image2D>& img_v, 
    const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& thrumu_v, 
    const std::vector< std::vector< const larcv::Pixel2D* > >& endpts_v, std::vector<int>& endpts_used  ) {

    PassOutput_t output;

    // I should probably use a chain-of-command pattern
    // here the input and output are specific to each submethod. there is no ambiguity on order of methods
    // maybe later...

    // make untagged image and cluster
    extractBaseClusters( passcfg, img_v, thrumu_v, endpts_v, output );

    // group clusters into secondary clusters
    findClusterLinks( passcfg, img_v, output );

    // we use the clusters to run astar
    analyzeClusters( passcfg, img_v, badch_v, thrumu_v, endpts_v, endpts_used, output );

    return output;
  }

  void StopMuCluster::extractBaseClusters( const StopMuClusterConfig::PassConfig_t& passcfg, const std::vector<larcv::Image2D>& img_v, 
    const std::vector<larcv::Image2D>& thrumu_v, const std::vector< std::vector< const larcv::Pixel2D* > >& endpts, StopMuCluster::PassOutput_t& output ) {

    // input checks
    if ( thrumu_v.size()!=img_v.size()) {
      throw std::runtime_error( "StopMuCluster::extractBaseClusters[Error] number of original and thrumu-tagged images are not the same.");
    }
    if ( img_v.size()!=m_config.pixel_thresholds.size() )
      throw std::runtime_error( "StopMuCluster::extractBaseClusters[Error] number of thresholds and images not the same.");

    // first mask out thrumu images
    output.m_masked_v.clear();
    for ( size_t iimg=0; iimg<img_v.size(); iimg++ ) {
      const larcv::Image2D& img    = img_v.at(iimg);
      const larcv::Image2D& thrumu = thrumu_v.at(iimg);

      // check that dimensions match
      if ( img.meta().rows()!=thrumu.meta().rows() || img.meta().cols()!=thrumu.meta().cols() ) {
        throw std::runtime_error("StopMuCluster::extractBaseClusters[Error] thrumu and orignal image dimensions do not match.");
      }

      larcv::Image2D masked( img.meta() );
      masked.paint(0.0);

      for (size_t r=0; r<img.meta().rows(); r++) {
        for (size_t c=0; c<img.meta().cols(); c++) {
          // skip if below threshold
          if ( img.pixel(r,c)<m_config.pixel_thresholds.at(iimg) ) 
            continue;
          // skip if tagged
          if ( thrumu.pixel(r,c)>0 )
            continue;
          masked.set_pixel(r,c,img.pixel(r,c));
        }
      }

      output.m_masked_v.emplace_back( std::move(masked) );


    }//end of image loop

    // unmask pixels around end points, so they can be included in the clusters
    for ( auto const& endpt : endpts ) {
      for (size_t p=0; p<output.m_masked_v.size(); p++ ) {
        larcv::Image2D& masked = output.m_masked_v.at(p);
        int row = (int)endpt.at(p)->Y();
        int col = (int)endpt.at(p)->X();
        for (int dr=-m_config.start_point_pixel_neighborhood; dr<m_config.start_point_pixel_neighborhood; dr++) {
          int r = row+dr;
          if ( r<0 || r>=(int)masked.meta().rows()) continue;
          for (int dc=-m_config.start_point_pixel_neighborhood; dc<m_config.start_point_pixel_neighborhood; dc++) { 
            int c = col+dc;
            if ( c<0 || c>=(int)masked.meta().cols() ) continue;
            masked.set_pixel(r,c,m_config.pixel_thresholds[p]+1);
          }
        }
      }
    }

    // use the masked image to form clusters
    output.m_untagged_clusters_v.clear();
    for (size_t p=0; p<output.m_masked_v.size(); p++) {
      untagged_cluster_info_t plane_cluster;
      plane_cluster.pixels = dbscan::extractPointsFromImage( output.m_masked_v.at(p), 0.5 );
      dbscan::DBSCANAlgo algo;
      plane_cluster.output = algo.scan( plane_cluster.pixels, m_config.dbscan_cluster_minpoints, m_config.dbscan_cluster_radius, false );
      for (size_t ic=0; ic<plane_cluster.output.clusters.size(); ic++) {
        dbscan::ClusterExtrema ex = dbscan::ClusterExtrema::FindClusterExtrema( ic, plane_cluster.output, plane_cluster.pixels );
        plane_cluster.extrema_v.emplace_back( std::move(ex) );
      }

      output.m_untagged_clusters_v.emplace_back( std::move(plane_cluster) );
    }

  }

  void StopMuCluster::findClusterLinks( const StopMuClusterConfig::PassConfig_t& passcfg, const std::vector<larcv::Image2D>& img_v, PassOutput_t& data ) {

    const size_t nplanes = img_v.size();
    for (size_t p=0; p<nplanes; p++) {
      const dbscan::dbClusters& clusters = data.m_untagged_clusters_v.at(p).output.clusters;
      for (size_t ic_a=0; ic_a<clusters.size(); ic_a++ ) {

        if ( (int)clusters.at(ic_a).size()<m_config.dbscan_cluster_minpoints ) continue;

        for (size_t ic_b=ic_a+1; ic_b<clusters.size(); ic_b++ ) {

          if ( (int)clusters.at(ic_b).size()<m_config.dbscan_cluster_minpoints ) continue;

          const dbscan::ClusterExtrema& ex_a = data.m_untagged_clusters_v.at(p).extrema_v.at(ic_a);
          const dbscan::ClusterExtrema& ex_b = data.m_untagged_clusters_v.at(p).extrema_v.at(ic_b);

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

          if ( min_dist>passcfg.alldir_max_link_dist && (min_dist <0 || min_dist > passcfg.max_link_distance ) ) continue;

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
            coscluster_a += dir_min[v]*dir_a[v];
            coscluster_b += -dir_min[v]*dir_b[v];
          }

          // cosine must be above some value for both
          if ( coscluster_b<passcfg.min_link_cosine || coscluster_a<passcfg.min_link_cosine )
            continue;

          // define the link
          data.m_untagged_clusters_v.at(p).makeLink( ic_a, ic_b, min_exa, min_exb, min_dist );
          // std::cout << "plane " << p << ": make link between " << ic_a << " and " << ic_b 
          //   << " ex(a)=" << min_exa << " ex(b)=" << min_exb
          //   << " a=(" << ex_a.extrema(min_exa)[0] << "," << ex_a.extrema(min_exa)[1] << ") "
          //   << " b=(" << ex_b.extrema(min_exb)[0] << "," << ex_b.extrema(min_exb)[1] << ") "            
          //   << " dist=" << min_dist << std::endl;
        }
      }
    }//end loop over planes

  }

  void StopMuCluster::analyzeClusters( const StopMuClusterConfig::PassConfig_t& passcfg, const std::vector<larcv::Image2D>& img_v, 
    const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& thrumu_v, const std::vector< std::vector< const larcv::Pixel2D* > >& endpts_v, 
    std::vector<int>& endpts_used, StopMuCluster::PassOutput_t& data  ) {

    // here we build stopmu-tracks using cluster groups found by previous pass sub-methods
    // for each space point, given as a pixel on all 3 planes, we
    //   (1) find the matching cluster group
    //   (2) look for 3D spacepoints on the cluster group
    //   (3) take the space point furthest from the start
    //   (4) and attempt to form a path from start to goal using AStar3DAlgo
    //   (5) we communicate if we are successful through the endpts_used vector

    // prep
    // build base clusters and links between them
    data.m_spacepoints.clear();
    data.m_paths.clear();
    data.m_path_goalreached.clear();
    data.m_clustergroups.clear();
    data.m_endpt_index.clear();

    if ( endpts_used.size()!=endpts_v.size() ) {
      endpts_used.resize( endpts_v.size(), 0 );
      // if the size is correct, we use that as a way to skip endpts for whatever reason
    }

    // setup A* config
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;    

    struct SPdata_t {
      int idx;
      float dist;
      bool operator()( SPdata_t a, SPdata_t b ) {
        if ( a.dist>b.dist ) return true;
        return false;
      };
    } myspdata;    

    for (size_t iendpt=0; iendpt<endpts_v.size(); iendpt++ ) {
      if ( endpts_used.at(iendpt)==1 ) continue;     

      auto const& endpt = endpts_v.at(iendpt);
      std::cout << "[Start Point idx=" << iendpt << "] " << std::endl;
      std::cout << " tick=" << img_v.front().meta().pos_y( endpt.front()->Y() ) << std::endl;
      std::cout << " cols=" << endpt[0]->X() << " " << endpt[1]->X() << " " << endpt[2]->X() << std::endl;

      // get 3D position of spacepoint
      float tick_start = img_v.at(0).meta().pos_y( endpt.at(0)->Y() );
      float x_start    = (tick_start-3200.0)*cm_per_tick;
      std::vector<int> wid_start(3,0);
      for ( size_t p=0; p<endpt.size(); p++ ) {
        wid_start[p] = img_v.at(p).meta().pos_x( endpt.at(p)->X() );
      }
      int cross_start = 0;
      double tri_start = 0.;
      std::vector<float> start_zy;
      larcv::UBWireTool::wireIntersection( wid_start, start_zy, tri_start, cross_start );
      std::vector<float> start_pos(3,0.0);
      start_pos[0] = x_start;
      start_pos[1] = start_zy[1];
      start_pos[2] = start_zy[0];

      // on each plane, we build cluster group connected to start points
      std::vector<int> starting_cluster(img_v.size(),-1);

      // get cluster connected to start point
      bool cluster_on_all_planes = true;
      for (size_t p=0; p<img_v.size(); p++) {
        std::vector<double> testpoint(2);
        testpoint[0] = endpt.at(p)->X();
        testpoint[1] = endpt.at(p)->Y();
        starting_cluster[p] = data.m_untagged_clusters_v.at(p).output.findMatchingCluster( testpoint, data.m_untagged_clusters_v.at(p).pixels, 2.0 );
        if ( starting_cluster[p]<0 ) {
          cluster_on_all_planes = false;
          break;
        }
      }

      if ( !cluster_on_all_planes ) {
        // found no compatible cluster on all three planes for this end point, we skip it
        // we fill an empty path for this end point though
        std::vector<AStar3DNode> empty;
        data.m_paths.emplace_back( std::move(empty) );
        data.m_path_goalreached.push_back(0);
        data.m_endpt_index.push_back( iendpt );
        continue;
      }

      // for each plane. build cluster group.
      ClusterGroup_t clustergroup;
      for ( int p=0; p<(int)img_v.size(); p++ ) {
        PlaneClusterGroup_t plgroup;
        std::vector<int> cluster_history;
        cluster_history.push_back( starting_cluster[p] );
        plgroup.group.insert(starting_cluster[p]);
        // recursive function
        getNextLinkedCluster( data, p, cluster_history, plgroup.group, plgroup.links );
        clustergroup.emplace_back( std::move(plgroup) );
      }

      // with cluster group made, we fill an image with pixels, which we will pass to A*
      std::vector< larcv::Image2D > clust_img_v;
      for ( size_t p=0; p<img_v.size(); p++ ) {
        larcv::Image2D clust_img( img_v.at(p).meta() );
        clust_img.paint(0.0);
        for ( auto &idx_cl : clustergroup.at(p).group ) {
          //std::cout << "plane " << p << " cluster group: " << idx_cl << std::endl;
          const dbscan::dbCluster& cluster = data.m_untagged_clusters_v.at(p).output.clusters.at(idx_cl);
          for (size_t ihit=0; ihit<cluster.size(); ihit++) {
            int hitidx = cluster.at(ihit);
            int col = data.m_untagged_clusters_v.at(p).pixels.at(hitidx)[0];
            int row = data.m_untagged_clusters_v.at(p).pixels.at(hitidx)[1];
            clust_img.set_pixel( row, col, img_v.at(p).pixel(row,col) );
            std::vector<int> pix(2);
            pix[0] = col;
            pix[1] = row;
            clustergroup.at(p).pixels.emplace_back( std::move(pix) );
          }
        }

        // we also use the links to fill pixels as well
        for ( auto const& plink : clustergroup.at(p).links ) {
          // we step through, filling in empty pixels where we can
          int nsteps = (*plink).dist/m_config.link_stepsize+1;
          float stepsize = (*plink).dist/float(nsteps);
          const std::vector<double>& start = data.m_untagged_clusters_v.at(p).extrema_v.at( (*plink).indices[0] ).extrema( (*plink).extrema[0] );
          const std::vector<double>& end   = data.m_untagged_clusters_v.at(p).extrema_v.at( (*plink).indices[1] ).extrema( (*plink).extrema[1] );
          double dir[2] = { (end[0]-start[0])/(*plink).dist, (end[1]-start[1])/(*plink).dist };
          for ( int istep=0; istep<nsteps; istep++) {
            int col = start[0] + stepsize*dir[0]*(istep+1);
            int row = start[1] + stepsize*dir[1]*(istep+1);
            for (int dr=-m_config.start_point_pixel_neighborhood; dr<=m_config.start_point_pixel_neighborhood; dr++) {
              int r = row+dr;
              if ( r<0 || r>=(int)clust_img.meta().rows() ) continue;
              for (int dc=-m_config.start_point_pixel_neighborhood; dc<=m_config.start_point_pixel_neighborhood; dc++) {
                int c = col+dc;
                if ( c<0 || c>=(int)clust_img.meta().cols() ) continue;
                if ( clust_img.pixel(r,c)<m_config.pixel_thresholds[p] ) {
                  clust_img.set_pixel(r,c,2.0*m_config.pixel_thresholds[p]);
                  std::vector<int> pix(2);
                  pix[0] = c;
                  pix[1] = r;
                  clustergroup.at(p).pixels.emplace_back(std::move(pix));
                }
              }
            }
          }
        }

        clust_img_v.emplace_back( std::move(clust_img) );
      }//end of loop over planes

      // we find spacepoints
      std::vector<BoundarySpacePoint> cluster_spacepoints = generateCluster2PlaneSpacepoints( passcfg, clustergroup, img_v, badch_v, thrumu_v, data );
      std::cout << "cluster space points: " << cluster_spacepoints.size() << std::endl;
      if ( cluster_spacepoints.size()==0 ) {
        // no space points. fill empty path
        std::vector<AStar3DNode> empty;
        data.m_paths.emplace_back( std::move(empty) );
        data.m_path_goalreached.push_back(0);
        data.m_endpt_index.push_back( iendpt );        
        continue;
      }

      // compress image for A* 3D
      std::vector< larcv::Image2D > clust_img_compressed_v;
      std::vector< larcv::Image2D > badch_compressed_v;
      float ds = m_config.astar_downsampling_factor;

      for ( size_t p=0; p<clust_img_v.size(); p++) {
        auto const& clust_img = clust_img_v.at(p);
        auto const& badch     = badch_v.at(p);
        larcv::Image2D compressed( clust_img.meta() );
        larcv::Image2D badch_compressed( clust_img.meta() );
        badch_compressed.paint(0);

        for (size_t r=0; r<clust_img.meta().rows(); r++) {
          for (size_t c=0; c<clust_img.meta().cols(); c++) {
            compressed.set_pixel(r,c, clust_img.pixel(r,c));
            badch_compressed.set_pixel(r,c,badch.pixel(r,c));
          }
        }
        compressed.compress( clust_img.meta().rows()/ds, clust_img.meta().cols()/ds );
        badch_compressed.compress( clust_img.meta().rows()/ds, clust_img.meta().cols()/ds );

        clust_img_compressed_v.emplace_back( std::move(compressed) );
        badch_compressed_v.emplace_back( std::move(badch_compressed) );
      }

      // sort through boundary space points by distance from start point

      std::vector< SPdata_t > spdata_v;
      for ( size_t idx=0; idx<cluster_spacepoints.size(); idx++) {
        const BoundarySpacePoint& sp = cluster_spacepoints.at(idx);
        // dist to start
        float dist = 0.;
        for (int v=0; v<3; v++) {
          float dv = sp.pos()[v]-start_pos[v];
          dist += dv*dv;
        }
        dist = sqrt(dist);
        SPdata_t spdata;
        spdata.idx = idx;
        spdata.dist = dist;
        spdata_v.emplace_back( std::move(spdata) );
      }

      // sort from smallest to largest
      std::sort( spdata_v.begin(), spdata_v.end(), myspdata );

      // we try to fit path from start to many of the end-points

      // first, get the start point
      int start_row    = endpt.at(0)->Y();
      float start_tick = img_v.at(0).meta().pos_y( start_row );
      std::vector<int> start_cols(3,0);
      std::vector<int> start_wids(3,0);
      for (size_t p=0; p<img_v.size(); p++) {
        start_cols[p] = endpt.at(p)->X();
        start_wids[p] = img_v.at(p).meta().pos_x( endpt.at(p)->X() );
      }

      // the spacepoints should be sorted in distance order, furthest to closest
      std::vector<AStar3DNode> path;
      bool goodtrack = false;
      int spidx = 0;
      for ( auto const& spdata : spdata_v ) {

        // Define Goal from space point
        int goal_row    = cluster_spacepoints.at(spdata.idx).at(0).row;
        float goal_tick = img_v.front().meta().pos_y( goal_row );
        std::vector<int> goal_cols(3,0);
        std::vector<int> goal_wids(3,0);
        for (size_t p=0; p<img_v.size(); p++) {
          goal_cols[p] = cluster_spacepoints.at(spdata.idx).at(p).col;
          goal_wids[p] = img_v.at(p).meta().pos_x( cluster_spacepoints.at(spdata.idx).at(p).col );
        }


        // std::cout << "start tick=" << start_tick << "; "
        //           << " start row=" << start_row  << "; "
        //           << " start col: " << endpt.at(0)->X() << " " << endpt.at(1)->X() << " " << endpt.at(2)->X() << "; "
        //           << std::endl;
        std::cout << "goal tick=" << goal_tick << "; goal row=" << goal_row << "; " 
                  << "goal cols: " << goal_cols[0] << " " << goal_cols[1] << " " << goal_cols[2] << std::endl;


        // first run linear tracker 
        goodtrack = false;
        PointInfoList linearpath = runLinearFitter( passcfg, img_v, badch_v, start_row, goal_row, start_cols, goal_cols, goodtrack );

        if ( goodtrack ) {
          // convert to vector<AStarNode> to make the same format as AStar
          std::vector<AStar3DNode> apath;
          for (int istep=(int)linearpath.size()-1; istep>=0; istep--) {
            AStar3DNode node;
            node.tyz = linearpath.at(istep).xyz;
            node.tyz[0] = node.tyz[0]/cm_per_tick+3200.0; // turn x into tick
            apath.emplace_back( std::move(node) );
          }
          std::swap( path, apath );
        }
        if ( (goodtrack && (linearpath.fractionGood()<0.9 || linearpath.fractionHasChargeOnMajorityOfPlanes()<0.9))
            || ( !goodtrack && (linearpath.fractionGood()>0.5 && linearpath.fractionHasChargeOnMajorityOfPlanes()>0.5) ) ) {
          // if not good, we then try the astar tracker
          bool goodastar = false;
          std::vector<AStar3DNode> apath = runAStar( passcfg, clust_img_v, badch_v, clust_img_compressed_v, badch_compressed_v, 
            start_row, goal_row, start_cols, goal_cols, goodastar );
          if ( goodastar ) {
            goodtrack = goodastar;
            std::swap(path,apath);
          }
        }

        // store spacepoints for drawing
        //for ( auto& sp : cluster_spacepoints ) 
        //m_spacepoints.emplace_back( std::move(sp) );
        data.m_spacepoints.emplace_back( std::move(cluster_spacepoints.at(spdata.idx)) );

        if ( goodtrack )
          break;
        spidx++;
        //if ( spidx>20 )
        //  break;
      }

      if ( !goodtrack )
        path.clear();
      data.m_paths.emplace_back( std::move(path) );

      if (goodtrack)
        data.m_path_goalreached.push_back( 1 );
      else
        data.m_path_goalreached.push_back( 0 );        
      data.m_endpt_index.push_back( iendpt );

      if ( goodtrack )
        endpts_used[iendpt] = 1;

      // store cluster object
      data.m_clustergroups.emplace_back( std::move(clustergroup) );

    }// end of endpoint loop

  }

  std::vector<BoundarySpacePoint> StopMuCluster::generateCluster3PlaneSpacepoints( const StopMuClusterConfig::PassConfig_t& passcfg, 
    const std::vector<larcv::Image2D>& img_v, const std::vector< std::set<int> >& cluster_groups, StopMuCluster::PassOutput_t& data ) {
    // get cluster group from 3 planes. collect the consistent 3 plane spacepoints

    // this involves a p*(p-1) combinatoric search...
    // we exclude by time quickly, before attempting intersection test

    // collect the different extrema points on the plane for easier processing
    std::cout << "cluster group 3-plane points" << std::endl;
    typedef std::vector< std::vector<float> > Pointlist_t;
    std::vector< Pointlist_t > planepts;
    for (size_t p=0; p<cluster_groups.size(); p++) {
      Pointlist_t list;
      std::cout << " plane " << p << ": ";
      for ( auto& cl_idx : cluster_groups.at(p)) {
        std::cout << cl_idx << " ";
        const dbscan::ClusterExtrema& ex = data.m_untagged_clusters_v.at(p).extrema_v.at(cl_idx);
        for (int i=0; i<(int)dbscan::ClusterExtrema::kNumExtrema; i++) {
          std::vector<float> pt(2);
          for (int v=0; v<2; v++)
            pt[v] = ex.extrema( (dbscan::ClusterExtrema::Extrema_t)i )[v];
          list.emplace_back( std::move(pt) );
        }
      }
      std::cout << std::endl;
      planepts.emplace_back( std::move(list) );
    }

    // do some gross combinatorics...and find spacepoints
    std::vector<BoundarySpacePoint> spacepoints;
    for (int i=0; i<(int)planepts.at(0).size(); i++) {

      const std::vector<float>& pt_i = planepts.at(0).at(i);

      for (int j=0; j<(int)planepts.at(1).size(); j++) {

        const std::vector<float>& pt_j = planepts.at(1).at(j);
        if ( fabs(pt_i[1]-pt_j[1])>passcfg.max_extrema_row_diff )
          continue;

        for (int k=0; k<(int)planepts.at(2).size(); k++) {

          const std::vector<float>& pt_k = planepts.at(2).at(k);
          if ( fabs(pt_i[1]-pt_k[1])>passcfg.max_extrema_row_diff || fabs(pt_j[1]-pt_k[1])>passcfg.max_extrema_row_diff )
            continue;

          // ok, all within some acceptable time
          std::vector<int> wids(3);
          wids[0] = img_v.at(0).meta().pos_x( pt_i[0] );
          wids[1] = img_v.at(1).meta().pos_x( pt_j[0] );
          wids[2] = img_v.at(2).meta().pos_x( pt_k[0] );

          double tri = 0;
          int crosses = 0;
          std::vector<float> intersection_zy;
          larcv::UBWireTool::wireIntersection( wids, intersection_zy, tri, crosses );

          // get tick
          float ave_row = (pt_i[1]+pt_j[1]+pt_k[1])/3.0;
          float tick = img_v.at(0).meta().pos_y(ave_row);

          std::cout << " tick=" << tick << " wids: " << wids[0] << " " << wids[1] << " " << wids[2] << " tri=" << tri << std::endl;

          if ( crosses==0 || tri>passcfg.max_extrema_triarea )
            continue;

          float xpos = (tick-3200.0)*0.5*::larutil::LArProperties::GetME()->DriftVelocity();

          // acceptable space point
          BoundarySpacePoint sp;
          sp.push_back( BoundaryEndPt( pt_i[0], pt_i[1] )); 
          sp.push_back( BoundaryEndPt( pt_j[0], pt_j[1] )); 
          sp.push_back( BoundaryEndPt( pt_k[0], pt_k[1] ));

          std::vector<float> pos(3,0.0);
          pos[0] = xpos;
          pos[1] = intersection_zy[1];
          pos[2] = intersection_zy[0];
          sp.pos( pos );

          spacepoints.emplace_back( std::move(sp) );
        }//end of k loop
      }//end of j loop
    }//end of i loop
    return spacepoints;
  }

  std::vector<BoundarySpacePoint> StopMuCluster::generateCluster2PlaneSpacepoints( const StopMuClusterConfig::PassConfig_t& passcfg, 
    const ClusterGroup_t& cluster_group, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
    const std::vector<larcv::Image2D>& thrumu_v, StopMuCluster::PassOutput_t& data ) {

    struct PlanePoint_t {
      std::vector<float> pt;
      int plane;
      PlanePoint_t( float x, float y, int p ) {
        pt.resize(2);
        pt[0] = x;
        pt[1] = y;
        plane = p;
      };
    };

    std::cout << "Collect Plane Points for cluster group: " << std::endl;
    std::vector<PlanePoint_t> planepts;
    for (size_t p=0; p<cluster_group.size(); p++) {
      std::cout << " plane " << p << ": ";
      for ( auto& cl_idx : cluster_group.at(p).group ) {
        std::cout << cl_idx << " ";
        const dbscan::ClusterExtrema& ex = data.m_untagged_clusters_v.at(p).extrema_v.at(cl_idx);
        for (int i=0; i<(int)dbscan::ClusterExtrema::kNumExtrema; i++) {
          std::vector<float> pt(2);
          for (int v=0; v<2; v++)
            pt[v] = ex.extrema( (dbscan::ClusterExtrema::Extrema_t)i )[v];
          PlanePoint_t ppt( pt[0], pt[1], p );
          planepts.emplace_back( std::move(ppt) );
        }
      }
      std::cout << std::endl;
    }
    std::cout << " total extrema points: " << planepts.size() << std::endl;

    std::vector<BoundarySpacePoint> spacepoints;
    for (int ipt=0; ipt<(int)planepts.size(); ipt++) {
      const PlanePoint_t& pt_i = planepts.at(ipt);
      for (int jpt=ipt+1; jpt<(int)planepts.size(); jpt++) {
        const PlanePoint_t& pt_j = planepts.at(jpt);
        if ( pt_i.plane==pt_j.plane ) continue; // not interest in same-plane matches

        if ( fabs(pt_i.pt[1]-pt_j.pt[1])>passcfg.max_extrema_row_diff ) continue; // too out of time

        int wire_i = img_v.at(pt_i.plane).meta().pos_x( pt_i.pt[0]);
        int wire_j = img_v.at(pt_j.plane).meta().pos_x( pt_j.pt[0]);

        // in-time
        int crosses = 0;
        std::vector<float> intersection_zy(2,0.0);
        int otherplane = -1;
        int otherwire  = -1;
        larcv::UBWireTool::getMissingWireAndPlane( pt_i.plane, wire_i, pt_j.plane, wire_j, otherplane, otherwire, intersection_zy, crosses );

        if ( crosses==0 ) {
          //std::cout << "  does not cross: p1=(" << pt_i.pt[0] << "," << pt_i.pt[1] << ")  p2=(" << pt_j.pt[0] << "," << pt_j.pt[1] << ")" << std::endl;
          continue;
        }

        // check if other point has charge or is in badch
        bool foundcharge = false;
        int averow = 0.5*( pt_i.pt[1] + pt_j.pt[1] );
        float tick = img_v.at(0).meta().pos_y(averow);
        float xpos = (tick-3200.0)*0.5*::larutil::LArProperties::GetME()->DriftVelocity();

        int col = img_v.at(otherplane).meta().col( otherwire );

        for (int dr=-m_config.start_point_pixel_neighborhood; dr<=m_config.start_point_pixel_neighborhood; dr++) {
          int r = averow+dr;
          if ( r<0 || r>=(int)img_v.at(otherplane).meta().rows() ) continue;
          for (int dc=-m_config.start_point_pixel_neighborhood; dc<=m_config.start_point_pixel_neighborhood; dc++) {
            int c = col+dc;
            if ( c<0 || c>=(int)img_v.at(otherplane).meta().cols() ) continue;
            if ( img_v.at(otherplane).pixel(r,c)>m_config.pixel_thresholds[otherplane] || badch_v.at(otherplane).pixel(r,c)>0 ) {
              foundcharge = true;
              break;
            }
          }
          if (foundcharge)
            break;
        }

        std::vector< std::vector<float> > pts(3);
        pts[ pt_i.plane ].resize(2);
        pts[ pt_i.plane ][0] = pt_i.pt[0];
        pts[ pt_i.plane ][1] = averow;
        pts[ pt_j.plane ].resize(2);
        pts[ pt_j.plane ][0] = pt_j.pt[0];
        pts[ pt_j.plane ][1] = averow;
        pts[ otherplane ].resize(2);
        pts[ otherplane ][0] = col;
        pts[ otherplane ][1] = averow;

        if ( !foundcharge ) continue;        

        // std::cout << "  candidate 2-plane sp: tick=" << tick << " averow=" << averow
        //   << " wids=(" << pts[0][0] << "," << pts[1][0] << "," << pts[2][0] << ") "
        //   << " foundcharge=" << foundcharge
        //   << std::endl;        

        // define space point
        BoundarySpacePoint sp;
        for (int p=0; p<3; p++) {
          sp.push_back( BoundaryEndPt( pts[p][1], pts[p][0] ) );
        }

        std::vector<float> pos(3,0.0);
        pos[0] = xpos;
        pos[1] = intersection_zy[1];
        pos[2] = intersection_zy[0];
        sp.pos( pos );

        if ( std::find( spacepoints.begin(), spacepoints.end(), sp )==spacepoints.end() )
          spacepoints.emplace_back( std::move(sp) );
      }//end of j loop
    }//end of i loop

    std::cout << "number of 2-plane spacepoints: " << spacepoints.size() << std::endl;
    return spacepoints;
  }

  void StopMuCluster::getNextLinkedCluster( StopMuCluster::PassOutput_t& data, const int& plane, std::vector<int>& cluster_history, 
    std::set<int>& clustergroup, std::vector< const ClusterLink_t* >& used_links ) {

    // get the links
    if ( cluster_history.size()==0)
      return;
    int current_cluster = cluster_history.back();
    const std::vector<ClusterLink_t>& cluster_link = data.m_untagged_clusters_v.at(plane).getLinks(current_cluster);
    for ( auto& link : cluster_link ) {
      // go to first cluster not already in the group
      if ( link.indices[0]==current_cluster && clustergroup.find( link.indices[1] )==clustergroup.end() ) {
        clustergroup.insert( link.indices[1] );
        cluster_history.push_back( link.indices[1] );
        used_links.push_back( &link );
        getNextLinkedCluster( data, plane, cluster_history, clustergroup, used_links );
      }
      else if ( link.indices[1]==current_cluster && clustergroup.find( link.indices[0] )==clustergroup.end() ) {
        clustergroup.insert( link.indices[0] );
        cluster_history.push_back( link.indices[0] );
        used_links.push_back( &link );        
        getNextLinkedCluster( data, plane, cluster_history, clustergroup, used_links );
      }
    }
    // no link
    cluster_history.pop_back();
    return;
  }

  BMTrackCluster3D StopMuCluster::makeBMTrackCluster3D( const std::vector<AStar3DNode>& path, 
    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<const larcv::Pixel2D*>& start_pt ) {

    const larcv::ImageMeta& meta = img_v.front().meta();
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick] 
    const int nplanes = img_v.size();

    // fill track3d data
    BMTrackCluster3D track3d;

    // Start Point Information
    track3d.start_type = (larlitecv::BoundaryEnd_t) int(start_pt.front()->Intensity());
    track3d.row_start  = start_pt.front()->Y();
    track3d.tick_start = img_v.front().meta().pos_y( track3d.row_start );
    track3d.start_wire.resize(nplanes,0);
    track3d.start3D.resize(nplanes,0);
    for (int i=0; i<nplanes; i++) {
      track3d.start3D[i] = path.back().tyz[i];
      track3d.start_wire[i] = meta.pos_x( start_pt[i]->X() );
    }   
    track3d.start3D[0] = (track3d.start3D[0]-3200)*cm_per_tick;

    // End Point Information
    track3d.end_type   = larlitecv::kUndefined;
    track3d.tick_end   = path.front().tyz.at(0);
    track3d.row_end    = meta.row( track3d.tick_end );
    track3d.end_wire.resize(nplanes,0);
    track3d.end3D.resize(nplanes,0);
    track3d.end3D[0] = (track3d.end3D[0]-3200)*cm_per_tick;    
    for (int i=1; i<nplanes; i++)
      track3d.end3D[i] = path.front().tyz[i];
    Double_t xyz_end[3];
    for (int i=0; i<nplanes; i++) 
      xyz_end[i] = track3d.end3D[i];    
    for (int i=0; i<nplanes; i++)
      track3d.end_wire[i] = (int)larutil::Geometry::GetME()->WireCoordinate(xyz_end,i);


    // Prepare Track2D objects and an empty image to track which pixels we've marked
    std::vector<larcv::Image2D> tagged_v;
    for (int p=0; p<nplanes; p++) {
      BMTrackCluster2D track2d;
      BoundaryEndPt start( track3d.row_start, meta.col( track3d.start_wire[p]), track3d.start_type );
      BoundaryEndPt end( track3d.row_end, meta.col( track3d.end_wire[p] ), larlitecv::kUndefined );
      track2d.start = start;
      track2d.end   = end;
      track2d.plane = p;
      track3d.plane_paths.emplace_back( std::move(track2d) );
      larcv::Image2D tagged( img_v.at(p).meta() );
      tagged.paint(0.0);
      tagged_v.emplace_back( std::move(tagged) );
    }

    float nbad_nodes = 0;
    float total_nodes = 0;
    int nnodes = (int)path.size();
    for ( int inode=nnodes-1; inode>=1; inode-- ) {

      const AStar3DNode& node      = path.at(inode);
      const AStar3DNode& next_node = path.at(inode-1);
      if ( node.badchnode )
        nbad_nodes+=1.0;
      total_nodes+=1.0;

      float dir3d[3];
      float step0[3];
      float dist = 0.;
      for (int i=0; i<3; i++) {
        dir3d[i] = next_node.tyz[i] - node.tyz[i];
        step0[i] = node.tyz[i];
      }
      dir3d[0] *= cm_per_tick;
      step0[0] = (step0[0]-3200.0)*cm_per_tick;
      for (int i=0; i<3; i++)
        dist += dir3d[i]*dir3d[i];
      dist = sqrt(dist);
      for (int i=0; i<3; i++)
        dir3d[i] /= dist;

      int nsteps = dist/m_config.link_stepsize+1;
      float stepsize = dist/float(nsteps);

      for (int istep=0; istep<=nsteps; istep++) {
        Double_t xyz[3];
        std::vector<double> pt(3,0.0);
        for (int i=0; i<3; i++) {
          xyz[i] = step0[i] + stepsize*istep*dir3d[i];
          pt[i] = xyz[i];
        }
        float tick = xyz[0]/cm_per_tick + 3200.0;
        if ( tick<=meta.min_y() || tick>=meta.max_y() ) continue;
        track3d.path3d.emplace_back( std::move(pt) );        
        int row = meta.row( tick );        
        std::vector<int> cols(3);
        for (int p=0; p<nplanes; p++) {
          cols[p] = meta.col( larutil::Geometry::GetME()->WireCoordinate(xyz,p) );

          for (int dr=-m_config.start_point_pixel_neighborhood; dr<=m_config.start_point_pixel_neighborhood; dr++) {
            int r = row+dr;
            if ( r<0 || r>=(int)meta.rows()) continue;
            for (int dc=-m_config.start_point_pixel_neighborhood; dc<=m_config.start_point_pixel_neighborhood; dc++) {
              int c = cols[p]+dc;
              if ( c<0 || c>=(int)meta.cols()) continue;
              // tag pixels that are (1) untagged && (2) above threshold or bad channels
              if ( tagged_v.at(p).pixel(r,c)==0 && (img_v.at(p).pixel(r,c)>m_config.pixel_thresholds.at(p) || badch_v.at(p).pixel(r,c)>0 ) ) {
                tagged_v.at(p).set_pixel(r,c,255);
                larcv::Pixel2D pix(c,r);
                track3d.plane_paths.at(p).pixelpath.emplace_back( std::move(pix) );
              }
            }
          }//end of row loop
        }//end of plane loop
      }//end of steps

    }//end of node loop

    return track3d;
  }

  PointInfoList StopMuCluster::runLinearFitter( const StopMuClusterConfig::PassConfig_t& passcfg, const std::vector<larcv::Image2D>& img_v, 
    const std::vector<larcv::Image2D>& badch_v, 
    const int start_row, const int goal_row, const std::vector<int>& start_cols, const std::vector<int>& goal_cols, bool& goodpath ) {

    // linear track
    Linear3DFitterConfig cfg;
    cfg.step_size = 3.0;
    Linear3DFitter lineartrack( cfg );
    PointInfoList path = lineartrack.findpath( img_v, badch_v, start_row, goal_row, start_cols, goal_cols );

    std::cout << "linear track result: goodfrac=" << path.fractionGood() << " majfrac=" << path.fractionHasChargeOnMajorityOfPlanes() << std::endl;


    // now we need to decide the quality of the fit
    goodpath = false;
    if ( path.fractionGood()>0.8 && path.fractionHasChargeOnMajorityOfPlanes()>0.8 )
      goodpath = true;

    return path;
  }

  std::vector<AStar3DNode> StopMuCluster::runAStar( const StopMuClusterConfig::PassConfig_t& passcfg, const std::vector<larcv::Image2D>& img_v, 
    const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& img_compressed_v, const std::vector<larcv::Image2D>& badch_compressed_v, 
    const int start_row, const int goal_row, const std::vector<int>& start_cols, const std::vector<int>& goal_cols, bool& goodpath ) {

    const larcv::ImageMeta& meta      = img_v.at(0).meta();
    const larcv::ImageMeta& meta_comp = img_compressed_v.at(0).meta();

    int start_row_compressed = meta_comp.row( meta.pos_y(start_row) );
    int goal_row_compressed  = meta_comp.row( meta.pos_y(goal_row) );
    std::vector<int> start_cols_compressed(3,0);
    std::vector<int> goals_cols_compressed(3,0);
    for (size_t p=0; p<img_v.size(); p++) {
      start_cols_compressed[p] = meta_comp.col( meta.pos_x(start_cols.at(p)) );
      goals_cols_compressed[p] = meta_comp.col( meta.pos_x(goal_cols.at(p)) );
    }
    if ( start_row_compressed<=0 ) start_row_compressed = 1;
    if ( start_row_compressed>=(int)img_compressed_v.front().meta().rows() ) start_row_compressed = (int)img_compressed_v.front().meta().rows() - 1;
    if ( goal_row_compressed<=0 ) goal_row_compressed = 1;
    if ( goal_row_compressed>=(int)img_compressed_v.front().meta().rows() )  goal_row_compressed  = (int)img_compressed_v.front().meta().rows() - 1;

    int goal_reached = 0;
    larlitecv::AStar3DAlgo algo( passcfg.astarcfg );
    algo.setVerbose(0);
        
    std::vector<AStar3DNode> path = algo.findpath( img_compressed_v, badch_compressed_v, badch_compressed_v, 
      start_row_compressed, goal_row_compressed, start_cols_compressed, goals_cols_compressed, goal_reached );

    std::cout << "astar result: goalreached=" << goal_reached << "path-length=" << path.size() << std::endl;

    // // FOR DEBUG
    //const larcv::Image2D& badch_compressed = badch_compressed_v.at(p);        
    // cv::Mat cvimg = larcv::as_mat_greyscale2bgr( compressed, 10, 500 );
    // for (int r=0; r<compressed.meta().rows(); r++) {
    //   for (int c=0; c<compressed.meta().cols(); c++) {
    //     if ( badch_compressed.pixel(r,c)>0) {
    //       cv::Vec3b& color = cvimg.at<cv::Vec3b>( cv::Point(c,r) );
    //       color[0] = 125;
    //       color[1] = 0;
    //       color[2] = 0;
    //     }
    //   }
    // }
    // cv::circle( cvimg, cv::Point(start_cols[p],start_row), 1, cv::Scalar(0,255,0), -1);
    // cv::circle( cvimg, cv::Point(goal_cols[p],goal_row), 1, cv::Scalar(0,0,255), -1);
    // std::stringstream ss;
    // ss << "test_cl" << iendpt << "_p" << p << ".jpg";
    // cv::imwrite( ss.str(), cvimg );

    if ( goal_reached==1 )
      goodpath = true;
    else
      goodpath = false;

    return path;
  }

  // ================================================================================================
  //  OPENCV FUNCTIONS
  // ================================================================================================

  void StopMuCluster::saveClusterImageOCV( std::string filename ) {
    // for visual evaluation, we dump out various information used/constructed by this class
#ifdef USE_OPENCV
    // std::vector<cv::Mat> cvimg_v = makeBaseClusterImageOCV();
    // for (size_t p=0; p<cvimg_v.size(); p++) {
    //   std::stringstream ss;
    //   ss << filename << "_p" << p << ".png";
    //   cv::imwrite( ss.str(), cvimg_v.at(p) );
    // }
#endif
  }

#ifdef USE_OPENCV
  std::vector<cv::Mat> StopMuCluster::makeBaseClusterImageOCV( const PassOutput_t& data, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v ) {

    // we draw an image, that highlights the clusters, links between them, and interesting space points
    // we drawn an image per plane

    TRandom rand(1);
    std::vector<cv::Mat> cvimgs_v;

    for ( size_t p=0; p<img_v.size(); p++ ) {
      const larcv::Image2D& img = img_v.at(p);
      const larcv::ImageMeta& meta = img.meta();
      const untagged_cluster_info_t& cluster_info = data.m_untagged_clusters_v.at(p);

      // first make a CV image we can have fun with
      cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img, m_config.pixel_thresholds[p], 100 );

      // color in the thrumu
      const larcv::Image2D& thrumu = thrumu_v.at(p);
      for (size_t r=0; r<thrumu.meta().rows(); r++) {
        for (size_t c=0; c<thrumu.meta().cols(); c++) {
          if ( thrumu.pixel(r,c)>0 ){
            cv::Vec3b& pix = cvimg.at<cv::Vec3b>( cv::Point(c,r) );
            pix[0] = 255;
            pix[1] = 0;
            pix[2] = 0;
          }
        }
      }

      // ok, we now label base clusters with colors!
      for (size_t icluster=0; icluster<cluster_info.output.clusters.size(); icluster++){
        if ( cluster_info.output.clusters.at(icluster).size()<m_config.dbscan_cluster_minpoints ) continue;
        const dbscan::dbCluster& cluster = cluster_info.output.clusters.at(icluster);
        // pick a color
        cv::Vec3b color;
        color[0] = (int)(rand.Uniform()*255);
        color[1] = (int)(rand.Uniform()*255);
        color[2] = (int)(rand.Uniform()*255);
        for (size_t ihit=0; ihit<cluster_info.output.clusters.at(icluster).size(); ihit++ ) {
          int hitidx = cluster_info.output.clusters.at(icluster).at(ihit);
          int x = cluster_info.pixels.at(hitidx)[0];
          int y = cluster_info.pixels.at(hitidx)[1];
          cvimg.at<cv::Vec3b>( cv::Point(x,y) ) = color;
        }
        // color in the extrema
        const dbscan::ClusterExtrema& ex = cluster_info.extrema_v.at(icluster);
        cv::circle(cvimg, cv::Point(ex.leftmost()[0],   ex.leftmost()[1]),   3, cv::Scalar(0,255,0),-1);
        cv::circle(cvimg, cv::Point(ex.topmost()[0],    ex.topmost()[1]),    3, cv::Scalar(255,255,0),-1);          
        cv::circle(cvimg, cv::Point(ex.rightmost()[0],  ex.rightmost()[1]),  3, cv::Scalar(0,255,255),-1);          
        cv::circle(cvimg, cv::Point(ex.bottommost()[0], ex.bottommost()[1]), 3, cv::Scalar(255,0,255),-1);

        std::stringstream slabel;
        slabel << icluster;
        cv::putText(cvimg,slabel.str(),cv::Point(ex.leftmost()[0],ex.leftmost()[1]),cv::FONT_HERSHEY_PLAIN, 1, cv::Scalar(255,255,255) );
      }//end of cluster loop

      for ( auto const& link : cluster_info.link_v ) {
        int icluster_a = link.indices[0];
        int icluster_b = link.indices[1];
        int x1 = cluster_info.extrema_v.at(icluster_a).extrema(link.extrema[0])[0];
        int y1 = cluster_info.extrema_v.at(icluster_a).extrema(link.extrema[0])[1];
        int x2 = cluster_info.extrema_v.at(icluster_b).extrema(link.extrema[1])[0];
        int y2 = cluster_info.extrema_v.at(icluster_b).extrema(link.extrema[1])[1];
        cv::line( cvimg, cv::Point(x1,y1), cv::Point(x2,y2), cv::Scalar(0,0,255), 2 );
      }

      /*
      for ( auto const& clust_img : m_cluster_images ) {
        const larcv::Image2D& clust_img_p = clust_img.at(p);
        for (size_t r=0; r<clust_img_p.meta().rows(); r++) {
          for (size_t c=0; c<clust_img_p.meta().cols(); c++) {
            if ( clust_img_p.pixel(r,c)>m_config.pixel_thresholds[p] ) {
              cv::Vec3b& pixcol = cvimg.at<cv::Vec3b>( cv::Point(c,r) );
              pixcol[0] = 255;
              pixcol[1] = 255;
              pixcol[2] = 255;
            }
          }
        }
      }
      */

      for ( auto const& sp : data.m_spacepoints ) {
        const BoundaryEndPt& endpt = sp.at(p);
        cv::circle(cvimg, cv::Point(endpt.col,endpt.row),   5, cv::Scalar(255,255,255),-1);
      }

      /*
      std::cout << "number of cluster groups: " << m_clustergroups.size() << std::endl;
      for ( auto const& clgroup : m_clustergroups ) {
        for ( auto const& pix : clgroup.at(p).pixels ) {
          cv::Vec3b& pixcol = cvimg.at<cv::Vec3b>( cv::Point(pix[0],pix[1]) );
          pixcol[0] = 80;
          pixcol[1] = 127;
          pixcol[2] = 255;
        }
      }
      */

      for ( size_t ipath=0; ipath<data.m_paths.size(); ipath++ ) {
        auto const& path = data.m_paths.at(ipath);
        for ( int inode=0; inode<(int)(path.size()-1); inode++ ) {

          float tick_start = path.at(inode).tyz.at(0);
          Double_t xyz_start[3] = { (tick_start-3200.0)*0.5*0.110, path.at(inode).tyz.at(1), path.at(inode).tyz.at(2) };
          float wid_start  = larutil::Geometry::GetME()->WireCoordinate( xyz_start, p );

          float tick_end = path.at(inode+1).tyz.at(0);
          Double_t xyz_end[3] = { (tick_start-3200.0)*0.5*0.110, path.at(inode+1).tyz.at(1), path.at(inode+1).tyz.at(2) };
          float wid_end  = larutil::Geometry::GetME()->WireCoordinate( xyz_end, p );

          if ( data.m_path_goalreached.at(ipath) )
            cv::line( cvimg, cv::Point( meta.col(wid_start), meta.row(tick_start)), cv::Point( meta.col(wid_end), meta.row(tick_end) ), cv::Scalar(0,255,0), 3 );
          else
            cv::line( cvimg, cv::Point( meta.col(wid_start), meta.row(tick_start)), cv::Point( meta.col(wid_end), meta.row(tick_end) ), cv::Scalar(0,100,0), 3 );
        }
      }

      cvimgs_v.emplace_back( std::move(cvimg) );
    }

    return cvimgs_v;
  }
#endif  

}