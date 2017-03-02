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

namespace larlitecv {

  StopMuCluster::StopMuCluster( const StopMuClusterConfig& cfg ) : m_config(cfg) {
    // Constructor
    setVerbosity(m_config.verbosity);
  }

  void StopMuCluster::extractBaseClusters(const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v, 
    const std::vector< std::vector< const larcv::Pixel2D* > >& endpts ) {

    // input checks
    if ( thrumu_v.size()!=img_v.size()) {
      throw std::runtime_error( "StopMuCluster::extractBaseClusters[Error] number of original and thrumu-tagged images are not the same.");
    }
    if ( img_v.size()!=m_config.pixel_thresholds.size() )
      throw std::runtime_error( "StopMuCluster::extractBaseClusters[Error] number of thresholds and images not the same.");

    // first mask out thrumu images
    m_masked_v.clear();
    m_img_v.clear();
    m_thrumu_v.clear();
    for ( size_t iimg=0; iimg<img_v.size(); iimg++ ) {
      const larcv::Image2D& img    = img_v.at(iimg);
      const larcv::Image2D& thrumu = thrumu_v.at(iimg);

      // we store pointers to these items for future reference
      m_img_v.push_back( &img );
      m_thrumu_v.push_back( &thrumu );

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

      m_masked_v.emplace_back( std::move(masked) );


    }//end of image loop

    // unmask pixels around end points, so they can be included in the clusters
    for ( auto const& endpt : endpts ) {
      for (size_t p=0; p<m_masked_v.size(); p++ ) {
        larcv::Image2D& masked = m_masked_v.at(p);
        int row = (int)endpt.at(p)->Y();
        int col = (int)endpt.at(p)->X();
        for (int dr=-m_config.start_point_pixel_neighborhood; dr<m_config.start_point_pixel_neighborhood; dr++) {
          int r = row+dr;
          if ( r<0 || r>=masked.meta().rows()) continue;
          for (int dc=-m_config.start_point_pixel_neighborhood; dc<m_config.start_point_pixel_neighborhood; dc++) { 
            int c = col+dc;
            if ( c<0 || c>=masked.meta().cols() ) continue;
            masked.set_pixel(r,c,m_config.pixel_thresholds[p]+1);
          }
        }
      }
    }

    // use the masked image to form clusters
    m_untagged_clusters_v.clear();
    for (size_t p=0; p<m_masked_v.size(); p++) {
      untagged_cluster_info_t plane_cluster;
      plane_cluster.pixels = dbscan::extractPointsFromImage( m_masked_v.at(p), 0.5 );
      dbscan::DBSCANAlgo algo;
      plane_cluster.output = algo.scan( plane_cluster.pixels, m_config.dbscan_cluster_minpoints, m_config.dbscan_cluster_radius, false );
      for (size_t ic=0; ic<plane_cluster.output.clusters.size(); ic++) {
        dbscan::ClusterExtrema ex = dbscan::ClusterExtrema::FindClusterExtrema( ic, plane_cluster.output, plane_cluster.pixels );
        plane_cluster.extrema_v.emplace_back( std::move(ex) );
      }

      m_untagged_clusters_v.emplace_back( std::move(plane_cluster) );
    }

  }

  void StopMuCluster::findClusterLinks() {
    const size_t nplanes = m_img_v.size();
    for (size_t p=0; p<nplanes; p++) {
      const dbscan::dbClusters& clusters = m_untagged_clusters_v.at(p).output.clusters;
      for (size_t ic_a=0; ic_a<clusters.size(); ic_a++ ) {

        if ( m_untagged_clusters_v.at(p).output.clusters.at(ic_a).size()<m_config.dbscan_cluster_minpoints ) continue;

        for (size_t ic_b=ic_a+1; ic_b<clusters.size(); ic_b++ ) {

          if ( m_untagged_clusters_v.at(p).output.clusters.at(ic_b).size()<m_config.dbscan_cluster_minpoints ) continue;

          const dbscan::ClusterExtrema& ex_a = m_untagged_clusters_v.at(p).extrema_v.at(ic_a);
          const dbscan::ClusterExtrema& ex_b = m_untagged_clusters_v.at(p).extrema_v.at(ic_b);

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
          if ( coscluster_b<m_config.min_link_cosine || coscluster_a<m_config.min_link_cosine )
            continue;

          // define the link
          m_untagged_clusters_v.at(p).makeLink( ic_a, ic_b, min_exa, min_exb, min_dist );
          // std::cout << "plane " << p << ": make link between " << ic_a << " and " << ic_b 
          //   << " ex(a)=" << min_exa << " ex(b)=" << min_exb
          //   << " a=(" << ex_a.extrema(min_exa)[0] << "," << ex_a.extrema(min_exa)[1] << ") "
          //   << " b=(" << ex_b.extrema(min_exb)[0] << "," << ex_b.extrema(min_exb)[1] << ") "            
          //   << " dist=" << min_dist << std::endl;
        }
      }
    }//end loop over planes

  }

  void StopMuCluster::findStopMuTrack( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
    const std::vector<larcv::Image2D>& thrumu_v, const std::vector< std::vector< const larcv::Pixel2D* > >& endpts_v  ) {
    // here we build stopmu-tracks

    // prep
    // build base clusters and links between them
    m_cluster_images.clear();
    m_spacepoints.clear();
    m_paths.clear();

    // setup A* config
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    larlitecv::AStar3DAlgoConfig astar_config;
    astar_config.astar_threshold.resize(3,0);
    astar_config.astar_threshold[0] = 10.0;
    astar_config.astar_threshold[1] = 10.0;
    astar_config.astar_threshold[2] = 10.0;
    astar_config.astar_neighborhood.resize(3,4);
    astar_config.astar_start_padding = 4;
    astar_config.astar_end_padding   = 4;
    astar_config.lattice_padding = 10;
    astar_config.accept_badch_nodes = true;
    astar_config.min_nplanes_w_hitpixel = 3;
    larlitecv::AStar3DAlgo algo( astar_config );

    struct SPdata_t {
      int idx;
      float dist;
      bool operator()( SPdata_t a, SPdata_t b ) {
        if ( a.dist<b.dist ) return true;
        return false;
      };
    } myspdata;    

    for ( auto const& endpt : endpts_v ) {

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
        starting_cluster[p] = m_untagged_clusters_v.at(p).output.findMatchingCluster( testpoint, m_untagged_clusters_v.at(p).pixels, 2.0 );
        if ( starting_cluster[p]<0 ) {
          cluster_on_all_planes = false;
          break;
        }
      }

      if ( !cluster_on_all_planes )
        continue;

      //std::cout << starting_cluster[0] << " " << starting_cluster[1] << " " << starting_cluster[2] << std::endl;
      //if ( starting_cluster[2]!=37 )
      //  continue;

      // for each plane. build cluster group.
      std::vector< std::set<int> > cluster_groups;
      std::vector< std::vector< const ClusterLink_t* > > cluster_links;
      //std::vector< 
      for ( int p=0; p<(int)img_v.size(); p++ ) {
        std::set<int> group;
        std::vector< const ClusterLink_t* > links;
        std::vector<int> cluster_history;
        cluster_history.push_back( starting_cluster[p] );
        group.insert(starting_cluster[p]);
        // recursive function
        getNextLinkedCluster( p, cluster_history, group, links );
        cluster_groups.emplace_back( std::move(group) );
        cluster_links.emplace_back( std::move(links) );
      }

      // with cluster group made, we fill an image with pixels, which we will pass to A*
      std::vector< larcv::Image2D > clust_img_v;
      for ( size_t p=0; p<img_v.size(); p++ ) {
        larcv::Image2D clust_img( img_v.at(p).meta() );
        clust_img.paint(0.0);
        for ( auto &idx_cl : cluster_groups.at(p) ) {
          //std::cout << "plane " << p << " cluster group: " << idx_cl << std::endl;
          const dbscan::dbCluster& cluster = m_untagged_clusters_v.at(p).output.clusters.at(idx_cl);
          for (size_t ihit=0; ihit<cluster.size(); ihit++) {
            int hitidx = cluster.at(ihit);
            int col = m_untagged_clusters_v.at(p).pixels.at(hitidx)[0];
            int row = m_untagged_clusters_v.at(p).pixels.at(hitidx)[1];
            clust_img.set_pixel( row, col, img_v.at(p).pixel(row,col) );
          }
        }

        // we also use the links to fill pixels as well
        for ( auto const& plink : cluster_links.at(p) ) {
          // we step through, filling in empty pixels where we can
          int nsteps = (*plink).dist/m_config.link_stepsize+1;
          float stepsize = (*plink).dist/float(nsteps);
          const std::vector<double>& start = m_untagged_clusters_v.at(p).extrema_v.at( (*plink).indices[0] ).extrema( (*plink).extrema[0] );
          const std::vector<double>& end   = m_untagged_clusters_v.at(p).extrema_v.at( (*plink).indices[1] ).extrema( (*plink).extrema[1] );
          float dir[2] = { (end[0]-start[0])/(*plink).dist, (end[1]-start[1])/(*plink).dist };
          for ( int istep=0; istep<nsteps; istep++) {
            int col = start[0] + stepsize*dir[0]*(istep+1);
            int row = start[1] + stepsize*dir[1]*(istep+1);
            if ( col<0 || col>=clust_img.meta().cols() ) continue;
            if ( row<0 || row>=clust_img.meta().rows() ) continue;
            if ( clust_img.pixel(row,col)<m_config.pixel_thresholds[p] )
              clust_img.set_pixel(row,col,m_config.pixel_thresholds[p]+1);
          }
        }

        clust_img_v.emplace_back( std::move(clust_img) );
      }

      // we find spacepoints
      //std::vector<BoundarySpacePoint> cluster_spacepoints = generateCluster3PlaneSpacepoints( cluster_groups );  
      std::vector<BoundarySpacePoint> cluster_spacepoints = generateCluster2PlaneSpacepoints( cluster_groups, img_v, badch_v, thrumu_v );
      std::cout << "cluster space points: " << cluster_spacepoints.size() << std::endl;
      if ( cluster_spacepoints.size()==0 )
        continue;

      for ( auto& sp : cluster_spacepoints ) 
        m_spacepoints.emplace_back( std::move(sp) );

      // compress image for A* 3D
      std::vector< larcv::Image2D > clust_img_compressed_v;
      std::vector< larcv::Image2D > badch_compressed_v;
      float ds = m_config.astar_downsampling_factor;
      for ( auto const& clust_img : clust_img_v ) {
        larcv::Image2D compressed( clust_img.meta() );
        larcv::Image2D badch_compressed( clust_img.meta() );
        badch_compressed.paint(0);

        for (size_t r=0; r<clust_img.meta().rows(); r++) {
          for (size_t c=0; c<clust_img.meta().cols(); c++) {
            compressed.set_pixel(r,c, clust_img.pixel(r,c));
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

      std::sort( spdata_v.begin(), spdata_v.end(), myspdata );
      for (auto &data : spdata_v ) {
        std::cout << "  idx=" << data.idx << " dist=" << data.dist << std::endl;
      }

      // translate start and stop
      int start_row = clust_img_compressed_v.at(0).meta().row( img_v.at(0).meta().pos_y( endpt.at(0)->Y() ) );
      int goal_row  = clust_img_compressed_v.at(0).meta().row( img_v.at(0).meta().pos_y( cluster_spacepoints.at(spdata_v.back().idx).at(0).row ) );
      std::vector<int> start_cols(3,0);
      std::vector<int> goal_cols(3,0);
      for (int p=0; p<3; p++) {
        const larcv::Image2D& compressed = clust_img_compressed_v.at(p);
        const larcv::Image2D& img        = clust_img_v.at(p);
        start_cols[p] = compressed.meta().col( img.meta().pos_x( endpt.at(p)->X() ) );
        goal_cols[p]  = compressed.meta().col( img.meta().pos_x( cluster_spacepoints.at(spdata_v.back().idx).at(p).col ) );
      }
      int goal_reached = 0;
      std::cout << "start row=" << start_row << " tick=" << img_v.at(0).meta().pos_y( endpt.at(0)->Y() )
                << " start col: " << start_cols[0] << " " << start_cols[1] << " " << start_cols[2] << std::endl;
      std::cout << "goal row=" << goal_row << " tick=" << img_v.at(0).meta().pos_y( cluster_spacepoints.at(spdata_v.back().idx).at(0).row )
                << " goal col: " << goal_cols[0] << " " << goal_cols[1] << " " << goal_cols[2] << std::endl;

      std::vector<AStar3DNode> path = algo.findpath( clust_img_compressed_v, badch_compressed_v, badch_compressed_v, 
        start_row, goal_row, start_cols, goal_cols, goal_reached );
      std::cout << "astar attempt. goal-reached=" << goal_reached << " pathsize=" << path.size() << std::endl;
      m_paths.emplace_back( std::move(path) );

      // store cluster images
      m_cluster_images.emplace_back( std::move(clust_img_v) );
    }// end of endpoint loop

  }

  std::vector<BoundarySpacePoint> StopMuCluster::generateCluster3PlaneSpacepoints( const std::vector< std::set<int> >& cluster_groups ) {
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
        const dbscan::ClusterExtrema& ex = m_untagged_clusters_v.at(p).extrema_v.at(cl_idx);
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
        if ( fabs(pt_i[1]-pt_j[1])>m_config.max_extrema_row_diff )
          continue;

        for (int k=0; k<(int)planepts.at(2).size(); k++) {

          const std::vector<float>& pt_k = planepts.at(2).at(k);
          if ( fabs(pt_i[1]-pt_k[1])>m_config.max_extrema_row_diff || fabs(pt_j[1]-pt_k[1])>m_config.max_extrema_row_diff )
            continue;

          // ok, all within some acceptable time
          std::vector<int> wids(3);
          wids[0] = m_img_v.at(0)->meta().pos_x( pt_i[0] );
          wids[1] = m_img_v.at(1)->meta().pos_x( pt_j[0] );
          wids[2] = m_img_v.at(2)->meta().pos_x( pt_k[0] );

          double tri = 0;
          int crosses = 0;
          std::vector<float> intersection_zy;
          larcv::UBWireTool::wireIntersection( wids, intersection_zy, tri, crosses );

          // get tick
          float ave_row = (pt_i[1]+pt_j[1]+pt_k[1])/3.0;
          float tick = m_img_v.at(0)->meta().pos_y(ave_row);

          std::cout << " tick=" << tick << " wids: " << wids[0] << " " << wids[1] << " " << wids[2] << " tri=" << tri << std::endl;

          if ( crosses==0 || tri>m_config.max_extrema_triarea )
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

  std::vector<BoundarySpacePoint> StopMuCluster::generateCluster2PlaneSpacepoints( const std::vector< std::set<int> >& cluster_groups, 
    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& thrumu_v ) {

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

    std::cout << "Collect Plane Points for cluter group: " << std::endl;
    std::vector<PlanePoint_t> planepts;
    for (size_t p=0; p<cluster_groups.size(); p++) {
      std::cout << " plane " << p << ": ";
      for ( auto& cl_idx : cluster_groups.at(p)) {
        std::cout << cl_idx << " ";
        const dbscan::ClusterExtrema& ex = m_untagged_clusters_v.at(p).extrema_v.at(cl_idx);
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

        if ( fabs(pt_i.pt[1]-pt_j.pt[1])>m_config.max_extrema_row_diff ) continue; // too out of time

        int wire_i = m_img_v.at(pt_i.plane)->meta().pos_x( pt_i.pt[0]);
        int wire_j = m_img_v.at(pt_j.plane)->meta().pos_x( pt_j.pt[0]);

        // in-time
        int crosses = 0;
        std::vector<float> intersection_zy;
        int otherplane;
        int otherwire;
        larcv::UBWireTool::getMissingWireAndPlane( pt_i.plane, wire_i, pt_j.plane, wire_j, otherplane, otherwire, intersection_zy, crosses );

        if ( crosses==0 ) {
          //std::cout << "  does not cross: p1=(" << pt_i.pt[0] << "," << pt_i.pt[1] << ")  p2=(" << pt_j.pt[0] << "," << pt_j.pt[1] << ")" << std::endl;
          continue;
        }

        // check if other point has charge or is in badch
        bool foundcharge = false;
        int averow = 0.5*( pt_i.pt[1] + pt_j.pt[1] );
        float tick = m_img_v.at(0)->meta().pos_y(averow);
        float xpos = (tick-3200.0)*0.5*::larutil::LArProperties::GetME()->DriftVelocity();

        int col = m_img_v.at(otherplane)->meta().col( otherwire );

        for (int dr=-m_config.start_point_pixel_neighborhood; dr<=m_config.start_point_pixel_neighborhood; dr++) {
          int r = averow+dr;
          if ( r<0 || r>=m_img_v.at(otherplane)->meta().rows() ) continue;
          for (int dc=-m_config.start_point_pixel_neighborhood; dc<=m_config.start_point_pixel_neighborhood; dc++) {
            int c = col+dc;
            if ( c<0 || c>=m_img_v.at(otherplane)->meta().cols() ) continue;
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

        std::cout << "  candidate 2-plane sp: tick=" << tick << " averow=" << averow
          << " wids=(" << pts[0][0] << "," << pts[1][0] << "," << pts[2][0] << ") "
          << " foundcharge=" << foundcharge
          << std::endl;        

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

        spacepoints.emplace_back( std::move(sp) );
      }//end of j loop
    }//end of i loop

    std::cout << "number of 2-plane spacepoints: " << spacepoints.size() << std::endl;
    return spacepoints;
  }

  void StopMuCluster::getNextLinkedCluster( const int& plane, std::vector<int>& cluster_history, std::set<int>& clustergroup, 
    std::vector< const ClusterLink_t* >& used_links ) {

    // get the links
    if ( cluster_history.size()==0)
      return;
    int current_cluster = cluster_history.back();
    const std::vector<ClusterLink_t>& cluster_link = m_untagged_clusters_v.at(plane).getLinks(current_cluster);
    for ( auto& link : cluster_link ) {
      // go to first cluster not already in the group
      if ( link.indices[0]==current_cluster && clustergroup.find( link.indices[1] )==clustergroup.end() ) {
        clustergroup.insert( link.indices[1] );
        cluster_history.push_back( link.indices[1] );
        used_links.push_back( &link );
        getNextLinkedCluster( plane, cluster_history, clustergroup, used_links );
      }
      else if ( link.indices[1]==current_cluster && clustergroup.find( link.indices[0] )==clustergroup.end() ) {
        clustergroup.insert( link.indices[0] );
        cluster_history.push_back( link.indices[0] );
        used_links.push_back( &link );        
        getNextLinkedCluster( plane, cluster_history, clustergroup, used_links );
      }
    }
    // no link
    cluster_history.pop_back();
    return;
  }

  // ================================================================================================
  //  OPENCV FUNCTIONS
  // ================================================================================================

  void StopMuCluster::saveClusterImageOCV( std::string filename ) {
    // for visual evaluation, we dump out various information used/constructed by this class
#ifdef USE_OPENCV
    std::vector<cv::Mat> cvimg_v = makeBaseClusterImageOCV();
    for (size_t p=0; p<cvimg_v.size(); p++) {
      std::stringstream ss;
      ss << filename << "_p" << p << ".png";
      cv::imwrite( ss.str(), cvimg_v.at(p) );
    }
#endif
  }

#ifdef USE_OPENCV
  std::vector<cv::Mat> StopMuCluster::makeBaseClusterImageOCV() {

    // we draw an image, that highlights the clusters, links between them, and interesting space points
    // we drawn an image per plane

    TRandom rand(1);
    std::vector<cv::Mat> cvimgs_v;

    for ( size_t p=0; p<m_img_v.size(); p++ ) {
      const larcv::Image2D& img = *(m_img_v.at(p));
      const untagged_cluster_info_t& cluster_info = m_untagged_clusters_v.at(p);

      // first make a CV image we can have fun with
      cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img, m_config.pixel_thresholds[p], 100 );

      // color in the thrumu
      const larcv::Image2D& thrumu = *(m_thrumu_v.at(p));
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

      std::cout << "number of spacepoints: " <<  m_spacepoints.size() << std::endl;
      for ( auto const& sp : m_spacepoints ) {
        const BoundaryEndPt& endpt = sp.at(p);
        cv::circle(cvimg, cv::Point(endpt.col,endpt.row),   5, cv::Scalar(255,255,255),-1);
      }

      cvimgs_v.emplace_back( std::move(cvimg) );
    }

    return cvimgs_v;
  }
#endif  

}