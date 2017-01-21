#include "FlashMuonTaggerAlgo.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"

// larcv
#include "UBWireTool/UBWireTool.h"

// larlite
#include "LArUtil/Geometry.h"

namespace larlitecv {
    
  bool FlashMuonTaggerAlgo::flashMatchTrackEnds( const std::vector< larlite::event_opflash* >& opflashsets, 
                                                 const std::vector<larcv::Image2D>& tpc_imgs,
                                                 const std::vector<larcv::Image2D>& badch_imgs,
                                                 std::vector< std::vector< BoundaryEndPt > >& trackendpts ) {
    // output
    // ------
    //  trackendpts: list of vector of end pts. each (inner vector) is point on each plane.
    const int nplanes = tpc_imgs.size();
    trackendpts.clear();
    
    if ( nplanes==0 )
      return false;
    
    // get a meta
    const larcv::ImageMeta& meta = tpc_imgs.at(0).meta();
    
    // loop over all the flash containers
    for ( auto& ptr_event_flash : opflashsets ) {
      // loop over flashes
      for ( auto& opflash : *ptr_event_flash ) {
        
        // determine time of flash
        float tick_target = 0;
        float flash_tick = 0;
        BoundaryEndPt::BoundaryEnd_t point_type;
        std::string modename;
        if ( fSearchMode==kAnode ) {
          flash_tick = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick + 15;
          tick_target = flash_tick;
          point_type = BoundaryEndPt::kAnode;
          modename = "anode";
        }
        else if ( fSearchMode==kCathode ) {
          flash_tick = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
          tick_target = flash_tick + fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick+180.0;
          point_type = BoundaryEndPt::kCathode;          
          modename = "cathode";
        }
        else if ( fSearchMode==kOutOfImage ) {
          flash_tick  = opflash.Time();
          tick_target = flash_tick; // dummy opflash gives first or last tick of image
          point_type = BoundaryEndPt::kImageEnd;          
          modename = "image end";
        }
        else {
          //std::cout << "[ERROR] wrong search mode" << std::endl;
          return false;
        }
        
        // check if the opflash time occurs within the image
        if ( tick_target<meta.min_y() || tick_target>=meta.max_y() )
          continue;

        // first find the weighted mean and total q     
        float qtot = 0;
        float z_weighted = 0.;
        for (int ipmt=0; ipmt<32; ipmt++) {
          z_weighted += opflash.PE( ipmt )*pmtpos[ipmt][2];
          qtot += opflash.PE( ipmt );
        }
        if ( qtot>0 ) {
          z_weighted /= qtot;
        }

        // to set the range, we find the first hit above threshold from the mean
        float min_dist_z = 1e9;
        float max_dist_z = 0;
        for (int ipmt=0; ipmt<32; ipmt++) {
          float pe = opflash.PE(ipmt);
          float dist = pmtpos[ipmt][2]-z_weighted;
          if ( pe>5.0 ) {
            if ( dist<0 && min_dist_z>dist ) min_dist_z = dist;
            else if ( dist>0 && max_dist_z<dist ) max_dist_z = dist;
          }
        }

        // extend by some factor (example 10%). This is to ensure acceptance of tracks.
        float zwidth = fabs(min_dist_z)+fabs(max_dist_z);
        float extension = zwidth*fConfig.flash_zrange_extension*0.5;

        std::vector<float> z_range = { z_weighted+min_dist_z-extension, z_weighted+max_dist_z+extension}; // be more intelligent later
        std::vector<float> y_range = { -120.0, 120.0 };
        if ( fSearchMode==kOutOfImage ) {
          // accept all
          z_range[0] = 0;
          z_range[1] = 1100;
        }
        else {
          // just go to ends in these cases
          if ( z_range[0]<55.0 ) z_range[0] = 0; 
          if ( z_range[1]>980.0 ) z_range[1] = 1100;
        }
        
        int row_target = meta.row( tick_target );
        
        if ( fConfig.verbosity<1 ) {
          std::cout << "============================================================================================" << std::endl;
          std::cout << "[Opflash search] " << modename << " mode" << std::endl;
          std::cout << "  opflash: "
                    << " flash_tick=" << flash_tick
                    << " tick_target=" << tick_target
                    << " row_target=" << row_target
                    << " qtot= " << qtot
                    << " z_range=[" << z_range[0] << "," << z_range[1] << "] "
                    << " w_range=[" << int(z_range[0]/0.3/meta.pixel_width()) << "," << int(z_range[1]/0.3/meta.pixel_width()) << "] "
                    << " drift_t=" << fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick+240.0 << " ticks"
                    << std::endl;
        }
        
        // we find track ends on all three planes
        std::vector<dbscan::dbPoints> hits;
        for (int p=0; p<nplanes; p++) {
          const larcv::Image2D& img = tpc_imgs.at(p);
          dbscan::dbPoints planehits;
          for (int drow=-10; drow<=10; drow++) {
            int r_current = row_target + drow;
            if ( r_current<0 || r_current>=(int)meta.rows() )continue;
            for (int c=0; c<(int)meta.cols(); c++) {
              if ( img.pixel( r_current, c )>fConfig.pixel_value_threshold.at(p) ) {
                std::vector< double > pt(2);
                pt.at(0) = c;
                pt.at(1) = r_current;
                planehits.emplace_back(pt);
              }
            }
          }
          hits.emplace_back( std::move(planehits) );
        }
        
        // we cluster hits consistent with the flash
        std::vector<dbscan::dbscanOutput> dbscan_output; 
        std::vector< std::vector<ClusterInfo_t> > cluster_info;
        for (int p=0; p<nplanes; p++) {
          dbscan::DBSCANAlgo dbalgo;
          dbscan::dbscanOutput plane_dbscan_output      = dbalgo.scan( hits.at(p), 5, 5.0, false, 0.0 );
          std::vector<ClusterInfo_t> plane_cluster_info = analyzeClusters( plane_dbscan_output, hits.at(p), tpc_imgs.at(p), row_target, p, 5 );
          for ( auto& info : plane_cluster_info ) {
	    if ( fConfig.verbosity<1)
  	      std::cout << " plane=" << p << " flash-cluster: (w,t)=(" << tpc_imgs.at(p).meta().pos_x(info.col) << "," << tpc_imgs.at(p).meta().pos_y(info.row) << ")"
            	<< " hit-rmax=" << info.hits_rmax_boundary;
              if ( info.hits_rmax_boundary ) 
              	std::cout << " (" << tpc_imgs.at(p).meta().pos_x(info.col_max) << "," << tpc_imgs.at(p).meta().pos_y(info.row_max) << ")";
              else
              	std::cout << " (,)";
              std::cout << " hit-rmin=" << info.hits_rmin_boundary;
              if ( info.hits_rmin_boundary ) 
              	std::cout << " (" << tpc_imgs.at(p).meta().pos_x(info.col_min) << "," << tpc_imgs.at(p).meta().pos_y(info.row_min) << ")";
              else
              	std::cout << " (,)";
              std::cout << std::endl;
          }
          dbscan_output.emplace_back( std::move(plane_dbscan_output) );
          cluster_info.emplace_back( std::move(plane_cluster_info) );
        }

	std::vector< std::vector<BoundaryEndPt> > endpts_v;
        std::vector< std::vector<ClusterInfo_t> > accepted_cluster_matches;        
	findPlaneMatchedClusters( cluster_info, tpc_imgs, fConfig.max_triarea, z_range, y_range, endpts_v, accepted_cluster_matches );

	// transfer to output container
        if ( fConfig.verbosity<2 )
          std::cout << "Generated End Points" << std::endl;
        for ( auto& intersection : endpts_v ) {
          std::cout << "  (" << tpc_imgs.at(0).meta().pos_x(intersection[0].w) 
            << "," << tpc_imgs.at(1).meta().pos_x(intersection[1].w)
            << "," << tpc_imgs.at(2).meta().pos_x(intersection[2].w) << ")" 
            << std::endl;
          for ( auto& plane_pt : intersection )
	    plane_pt.type = point_type;
          trackendpts.emplace_back( std::move( intersection ) );
        }
        /*
        // now we loop through the clusters, and form track points on each plane
        // so that we define 'end points' on each plane        
        std::vector< std::vector<BoundaryEndPt> > endpts;        
        for ( int p=0; p<nplanes; p++ ) {
          const dbscan::dbscanOutput& plane_cluster_info = cluster_info.at(p);
          
          // for (int ic=0; ic<(int)plane_cluster_info.clusters.size(); ic++) {
          //   BoundaryEndPt endpt;
          //   bool foundend = findClusterEnds( plane_cluster_info.clusters.at(ic), hits.at(p), ic, row_target, p, meta, endpt );
          //   if ( foundend ) endpts[p]->emplace_back( endpt );
          // }
          
          std::vector<BoundaryEndPt> plane_endpts;
          defineClusterPositions( plane_cluster_info, hits.at(p), row_target, p, point_type, tpc_imgs.at(p), plane_endpts );
          endpts.emplace_back( std::move(plane_endpts) );
          if (fConfig.verbosity<1 )
            std::cout << "  flashmuontaggeralgo: plane " << p << " clusters: " << plane_cluster_info.clusters.size() << " endpoints=" << endpts[p].size() << std::endl;
        }

        // generate 3 plane intersections from clusters
        std::vector< std::vector<BoundaryEndPt> > intersections3plane = generate3PlaneIntersections( endpts, fConfig.max_triarea, tpc_imgs, z_range, y_range );

        // transfer to output container
        if ( fConfig.verbosity<2 )
          std::cout << "Generated End Points" << std::endl;
        for ( auto& intersection : intersections3plane ) {
          std::cout << "  (" << tpc_imgs.at(0).meta().pos_x(intersection[0].w) 
            << "," << tpc_imgs.at(1).meta().pos_x(intersection[1].w)
            << "," << tpc_imgs.at(2).meta().pos_x(intersection[2].w) << ")" 
            << std::endl;
          trackendpts.emplace_back( std::move( intersection ) );
        }
        */        

        /*
        // now that we have end points for each cluster, we ID position of charge deposition in (Y,Z) and check if consistent with the flash
        std::vector< std::vector<int> > wirelists(3);
        std::vector< std::map<int,int> > wirelist_endpt_idx(3); // map from wid to endpt index
        std::vector< std::vector<float> > valid_range(2);
        valid_range[0] = z_range;
        valid_range[1] = y_range;
        for (int ip=0; ip<3; ip++) {
          for (int ipt=0; ipt<(int)endpts[ip].size(); ipt++) {
            const BoundaryEndPt& endpt = endpts[ip].at(ipt);
            int wid = meta.pos_x(endpt.w);
            wirelists.at(ip).push_back( wid ); 
            wirelist_endpt_idx.at(ip).insert( std::pair<int,int>(wid,ipt) );
          }
        }
        
        // get intersections
        std::vector< std::vector<int> > intersections3plane;
        std::vector< std::vector<float> > vertex3plane;
        std::vector<float> areas3plane;
        std::vector< std::vector<int> > intersections2plane;
        std::vector< std::vector<float> > vertex2plane;
        
        ::larcv::UBWireTool::findWireIntersections( wirelists, valid_range, intersections3plane, vertex3plane, areas3plane, intersections2plane, vertex2plane );

        if ( fConfig.verbosity<1 ) {
          std::cout << " RAW INTERSECTIONS" << std::endl;
          std::cout << " 2 plane intersections: " << std::endl;
          for (int ii=0; ii<(int)intersections2plane.size(); ii++) {
            std::cout << "   (" << intersections2plane.at(ii).at(0) << ","
                    << intersections2plane.at(ii).at(1) << ","
                    << intersections2plane.at(ii).at(2) << ")"
                    << " vertex (z,y)=(" << vertex2plane.at(ii).at(0) << ", " << vertex2plane.at(ii).at(1) << ")" << std::endl;
          }
          
          std::cout << " 3 plane intersections:" << std::endl;
          for (int ii=0 ; ii<(int)intersections3plane.size(); ii++) {
            std::cout << "   (" << intersections3plane.at(ii).at(0) << ","
                    << intersections3plane.at(ii).at(1) << ","
                    << intersections3plane.at(ii).at(2) << ") "
                    << " vertex (z,y)=(" << vertex3plane.at(ii).at(0) << ", " << vertex3plane.at(ii).at(1) << ")"
                    << "tri-area-score=" << areas3plane.at(ii) << std::endl;
          }
        }

        // now we go through and find our best match(es)
        int max_flash_combos = 3;
        int ncombos = 0;
        std::vector< int > final_combos_idx;

        // ========================================================================
        // SORT 3 PLANE COMBINATIONS
        // sort indices by area
        struct p3matches_t {
          int idx;
          float score;
        };
        struct mycompare_t {
          bool operator() (p3matches_t l, p3matches_t r) {
            if ( l.score<r.score ) return true;
            return false;
          };
        } mycompare;

        std::vector< p3matches_t > sortable3plane;
        for (int ii=0; ii<(int)intersections3plane.size(); ii++) {
          p3matches_t combo;
          combo.idx = ii;
          combo.score = areas3plane.at(ii);
          sortable3plane.emplace_back( combo );
        }
        std::sort( sortable3plane.begin(), sortable3plane.end(), mycompare );

        // go through combos and pick out matches
        std::set<int> used_wid[3]; // used to mark wid's that are used up
        for (int ii=0; ii<(int)sortable3plane.size(); ii++) {
          const p3matches_t& combo = sortable3plane.at(ii);
          const std::vector<int>& wire_combo = intersections3plane.at( combo.idx );
          bool valid_combo = true;
//        for (int iw=0; iw<wire_combo.size(); iw++) {
//          int wid = wire_combo.at(iw);
//          if ( used_wid[iw].find( wid )!=used_wid[iw].end() ) {
//            // wire already used
//            valid_combo = false;
//            break;
//          }
//        }
          if ( !valid_combo || combo.score>fConfig.max_triarea )
            continue; // do not accept
          
          if ( ncombos>=max_flash_combos && combo.score>fConfig.max_triarea_tight )
            continue; // only accept more than the max unless its a really good match
          
          // valid combo. now invalid these wires
          for (int iw=0; iw<(int)wire_combo.size(); iw++) {
            int wid = wire_combo.at(iw);
            used_wid[iw].insert( wid );
          }
          
          final_combos_idx.push_back( combo.idx ); // copy combo
          ncombos++;
        }

        // =================================================================================================
        // SORT 2 PLANE COMBINATIONS: add in flash combinations where missing end point is in badch regions
        for ( int ii=0; ii<(int)intersections2plane.size(); ii++) {
          // find the missing one
          int missing_plane = -1;
          for (int ip=0; ip<3; ip++) {
            if ( intersections2plane.at(ii).at(ip)==-1 ) {
              missing_plane = ip;
              break;
            }
          }
          // get two plane intersection vertex
          std::vector<double> vert(3,0.0);
          vert[1] = vertex2plane.at(ii).at(1);
          vert[2] = vertex2plane.at(ii).at(0);
          float wireco = ::larutil::Geometry::GetME()->WireCoordinate( vert, missing_plane );
          // is it near a missing wire?
          int wirecol = (int)wireco/meta.pixel_width();
          bool missing_in_badch = false;
          for (int n=-3;n<=3;n++) {
            int col_ = wirecol+n;
            if ( col_<0 || col_>=(int)meta.cols() ) continue;
            if ( badch_imgs.at(missing_plane).pixel(0,col_)>0 ){
              missing_in_badch = true;
              break;
            }
          }
          if ( missing_in_badch ) {
            
            // add to 3 plane intersections
            std::vector<int> intersect3(3);
            for (int ip=0; ip<3; ip++)
              intersect3[ip] = intersections2plane.at(ii).at(ip);
            intersect3[missing_plane] = (int)wireco;
            int min_totdiff = 10000;
            for (int jj=0; jj<(int)intersections3plane.size(); jj++) {
              int totdiff = 0;
              for (int ip=0; ip<3; ip++)
                totdiff += abs( intersections3plane.at(jj).at(ip)-intersect3.at(ip) );
              if (totdiff<min_totdiff)
                min_totdiff = totdiff;
            }
            if ( min_totdiff>=5 ) {
              //std::cout << "  adding 2-plane intersection+badch wire into 3 plane intersections: "
              //<< "(" << intersect3[0] << "," << intersect3[1] << "," << intersect3[2] << ")" << std::endl;
              std::vector<float> fvertzy(2,0.0);
              fvertzy[0] = vertex2plane.at(ii).at(0);
              fvertzy[1] = vertex2plane.at(ii).at(1);
              int idx = intersections3plane.size();
              intersections3plane.emplace_back( std::move(intersect3) );
              vertex3plane.emplace_back( std::move(fvertzy) );
              areas3plane.push_back(0.0);
              final_combos_idx.push_back(idx);

              // add wire map entry if needed
              int wid = meta.pixel_width()*wirecol;
              if ( wirelist_endpt_idx.at(missing_plane).find( wid )==wirelist_endpt_idx.at(missing_plane).end() ) {
                wirelists.at(missing_plane).push_back( wid );
                
                BoundaryEndPt endpt( row_target, (int)wirecol);
                int idx_endpt = endpts[missing_plane].size();
                endpts[missing_plane].emplace_back( std::move(endpt) );
                wirelist_endpt_idx.at(missing_plane).insert( std::pair<int,int>(wid,idx_endpt) );
              }
            }
          }
        }
        

        // =================================================================================================
        // FINAL 3 PLANE COMBOS. FILL TRACK END POINTS CONTAINER WITH RESULTS
        if (fConfig.verbosity<2 )
          std::cout << "Final 3-plane combos" << std::endl;
        for (int ii=0 ; ii<(int)final_combos_idx.size(); ii++) {
          int idx = final_combos_idx.at(ii);
          std::vector<int> endptidx(3,-1);
          for (int ip=0; ip<nplanes; ip++) {
            int wid = intersections3plane.at(idx).at(ip);
            auto it_endpt = wirelist_endpt_idx.at(ip).find( wid );
            if ( it_endpt!=wirelist_endpt_idx.at(ip).end() ) {
              endptidx[ip] = it_endpt->second;
            }
          }
          if (fConfig.verbosity<2 ) {
            std::cout << "   (" << intersections3plane.at(idx).at(0) << ","
                      << intersections3plane.at(idx).at(1) << ","
                      << intersections3plane.at(idx).at(2) << ") "
                      << " endpt=(" << endptidx[0] << "," << endptidx[1] << "," << endptidx[2] << ") "
                      << " vertex (z,y)=(" << vertex3plane.at(idx).at(0) << ", " << vertex3plane.at(idx).at(1) << ")"
                      << " score=" << areas3plane.at(idx) << std::endl;
          }

          std::vector< BoundaryEndPt > endpt_v;
          for (int ip=0; ip<nplanes; ip++) {
            int eidx = endptidx.at(ip); // end point index
            larlitecv::BoundaryEndPt endpt = endpts[ip].at(eidx);
            if ( fSearchMode==kAnode )
              endpt.type = larlitecv::BoundaryEndPt::kAnode;
            else if ( fSearchMode==kCathode )
              endpt.type = larlitecv::BoundaryEndPt::kCathode;
            else if ( fSearchMode==kOutOfImage ) // makes no sense here
              endpt.type = larlitecv::BoundaryEndPt::kImageEnd;
            endpt_v.emplace_back( std::move( endpt ) );
          }//end of loop over planes
          // place into output container

          trackendpts.emplace_back( std::move( endpt_v ) );
        }
        */
        
      }//end of opflashes loop
    }//end of opflashsets

    return true;
  }
    
  void FlashMuonTaggerAlgo::loadGeoInfo() {

    TFile fGeoFile( Form("%s/app/PMTWeights/dat/geoinfo.root",getenv("LARCV_BASEDIR")), "OPEN" );

    // Get the PMT Info
    fNPMTs = 32;
    TTree* fPMTTree  = (TTree*)fGeoFile.Get( "imagedivider/pmtInfo" );
    int femch;
    float pos[3];
    fPMTTree->SetBranchAddress( "femch", &femch );
    fPMTTree->SetBranchAddress( "pos", pos );
    for (int n=0; n<fNPMTs; n++) {
      fPMTTree->GetEntry(n);
      for (int i=0; i<3; i++) {
        pmtpos[femch][i] = pos[i];
      }
      //std::cout << "[POS " << femch << "] " << " (" << pmtpos[femch][0] << "," << pmtpos[femch][1] << "," << pmtpos[femch][2] << ")" << std::endl;
    }
    
    fGeoFile.Close();

  }
  
  bool FlashMuonTaggerAlgo::findClusterEnds( const dbscan::dbscanOutput& clout, const dbscan::dbPoints& winpoints, 
                                             const int clusterid, const int row_target, const int plane, 
                                             const larcv::ImageMeta& meta, BoundaryEndPt& endpt   ) {
                                             
    
    // find extremities of cluster (in time)
    int tmax = -1;
    int tmin = -1;
    int wmax = -1;
    int wmin = -1;

    for ( int ichit=0; ichit<(int)clout.clusters.at(clusterid).size(); ichit++ ) {
      int hitidx = clout.clusters.at(clusterid).at(ichit);
      int x_ = (int)winpoints.at(hitidx).at(0)+0.1;
      int y_ = (int)winpoints.at(hitidx).at(1)+0.1;

      if ( tmax==-1 || y_>tmax ) { tmax = y_; wmax = x_; };
      if ( tmin==-1 || y_<tmin ) { tmin = y_; wmin = x_; };
    }
    if ( fConfig.verbosity<1 ) {
      std::cout << "  end points: max (r,c)=(" << tmax << ", " << wmax << ")"
                << "  tw=(" << meta.pos_y(tmax) << "," << meta.pos_x(wmax) << ")"
                << "  min (r,c)=(" << tmin << "," << wmin << ")" 
                << "  tw=(" << meta.pos_y(tmin) << "," << meta.pos_x(wmin) << ")"
                << "  versus: query_row=" << row_target
                << std::endl;
    }

    // is this a marked flash-end?
    // if extrema matches the annode flash hypothesis row. mark as interesting (Score 200)
    bool success = false;
    bool tmin_isend = false;
    bool tmax_isend = false;
    if ( abs(tmin-row_target)<=fConfig.endpoint_time_neighborhood.at(plane) ) {
      //if ( fConfig.verbosity<1 ) std::cout << "  MATCHED MIN END" << std::endl;
      tmin_isend = true;
    }
    if ( abs(tmax-row_target)<=fConfig.endpoint_time_neighborhood.at(plane) ) { 
      //if ( fConfig.verbosity<1 ) std::cout << "  MATCHED MAX_END" << std::endl;
      tmax_isend = true;
    }

    if ( (tmin_isend && !tmax_isend) || (tmin_isend && tmax_isend && abs(tmin-row_target)<abs(tmax-row_target)) ) {
      success = true;
      endpt.t = tmin;
      endpt.w = wmin;
      if ( fSearchMode==kAnode ) endpt.type = BoundaryEndPt::kAnode;
      else if ( fSearchMode==kCathode ) endpt.type = BoundaryEndPt::kCathode;
      else if ( fSearchMode==kOutOfImage ) endpt.type = BoundaryEndPt::kImageEnd;
    }
    else if ( (tmax_isend && !tmin_isend) || (tmax_isend && tmin_isend && abs(tmax-row_target)<abs(tmin-row_target)) ) {
      success = true;
      endpt.t = tmax;
      endpt.w = wmax;
      if ( fSearchMode==kAnode ) endpt.type = BoundaryEndPt::kAnode;
      else if ( fSearchMode==kCathode ) endpt.type = BoundaryEndPt::kCathode;
      else if ( fSearchMode==kOutOfImage ) endpt.type = BoundaryEndPt::kImageEnd;
    }      
        
    return success;
  }

  void FlashMuonTaggerAlgo::defineClusterPositions( const dbscan::dbscanOutput& clout, const dbscan::dbPoints& winpoints, const int row_target, 
    const int plane, const int point_type, const larcv::Image2D& img, std::vector<BoundaryEndPt>& endpt_v ) {
    // we generate a boundaryendpt for every cluster
    // we check the cluster size, is it big enough? 
    // we have a target row.  we scan across a cluster pixels and find the largest charge pixel on the row.
    // we use the above pixel to define the point
    // we can later condition on 3D-match points to recluster and find track ends
    int nclusters = clout.clusters.size();
    const larcv::ImageMeta& meta = img.meta();
    for ( int icluster=0; icluster<nclusters; icluster++ ) {
      const dbscan::dbCluster& cluster = clout.clusters.at(icluster);
      int nhits = cluster.size();
      if ( nhits < fConfig.flash_pixelcluster_minsize ) continue;
      // passes, find the pixel on row=row_target with the most charge
      float max_q = -1;
      int max_hitidx = -1;
      int max_col = -1;
      int min_row_gap = -1;
      int max_row = -1;
      for ( int ihit=0; ihit<nhits; ihit++ ) {
        int hitidx = cluster.at(ihit);
        const std::vector<double>& pixel = winpoints.at(hitidx);
        int row_gap = abs(row_target-(int)pixel[1]);
        float q = img.pixel( (int)pixel[1], (int)pixel[0] );        
        if ( min_row_gap<0 || min_row_gap>row_gap ) {
          // closest to row_target, set the min values
          max_q = q;          
          max_hitidx = hitidx;          
          max_col = (int)pixel[0];
          max_row = (int)pixel[1];
        }
        else if ( min_row_gap==row_gap ) {
          // same gap as we've seen before, so we pick which pixel has the highest charge
          // inspect pixels along the row target
          if ( max_q<0 || q>max_q ) {
            max_q = q;
            max_hitidx = hitidx;
            max_col = (int)pixel[0];
            max_row = (int)pixel[1];
          }
        }
      }
      // define the boundary point here
      if ( max_col>=0 ) {
        BoundaryEndPt endpt( max_row, max_col, (BoundaryEndPt::BoundaryEnd_t)point_type );
        endpt_v.emplace_back( std::move(endpt) );
      }
    }
  }

  std::vector< std::vector<BoundaryEndPt> > FlashMuonTaggerAlgo::generate3PlaneIntersections( const std::vector< std::vector<BoundaryEndPt> >&endptset_v, 
    const float max_triarea, const std::vector<larcv::Image2D>& img_v, const std::vector<float>& z_range, const std::vector<float>& y_range ) {

     if ( fConfig.verbosity<1){
        std::cout << "  generate 3-plane intersections within "
                << " z-range=[" << z_range[0] << "," << z_range[1] << "]" 
                << " y-range=[" << y_range[0] << "," << y_range[1] << "]" 
                << std::endl;           
     }

    // we try all combinations of 3 planes.  factorial expansion...
    std::vector<int> numpts(endptset_v.size());
    for (int p=0; p<(int)endptset_v.size(); p++) {
      numpts[p] = (int)endptset_v.at(p).size();
    }
    // dont feel like writing something that generates combinatorics, going to assume 3 planes
    std::vector< std::vector<int> > combos;
    for (int p1=0; p1<numpts[0]; p1++) {
      for (int p2=0; p2<numpts[1]; p2++) {
        for ( int p3=0; p3<numpts[2]; p3++ ) {
          std::vector<int> combo(3);
          combo[0] = p1;
          combo[1] = p2;
          combo[2] = p3;
          combos.emplace_back( std::move(combo) );
        }
      }
    }

    std::vector< std::vector<BoundaryEndPt> > output;
    for ( auto const& combo : combos ) {
      // get the wire ids for each combo
      std::vector<int> wid(endptset_v.size());
      for ( int p=0; p<3; p++) {
        wid[p] = img_v.at(p).meta().pos_x( endptset_v.at(p).at(combo[p]).w );
      }
      int crosses = 0;
      std::vector<float> poszy;
      double triangle_area;
      larcv::UBWireTool::wireIntersection( wid, poszy, triangle_area, crosses );
      if ( fConfig.verbosity<1)
        std::cout << "  candidate combo (" << wid[0] << "," << wid[1] << "," << wid[2] << ") "
                << " tri-area=" << triangle_area 
                << " poszy=(" << poszy[0] << "," << poszy[1] << ")"
                << " crosses=" << crosses 
                << std::endl;
      if ( crosses==1 && triangle_area < max_triarea 
            && z_range[0]<=poszy[0] && poszy[0]<=z_range[1] && y_range[0]<=poszy[1] && poszy[1]<=y_range[1] ) {
        std::vector<BoundaryEndPt> endpt_v;
        for (int p=0; p<3; p++) {
          endpt_v.push_back( endptset_v.at(p).at(combo[p]) );
        }
        output.emplace_back( std::move(endpt_v) );
      }
    }
    return output;
  }

  std::vector<FlashMuonTaggerAlgo::ClusterInfo_t> FlashMuonTaggerAlgo::analyzeClusters( const dbscan::dbscanOutput& dbscan_output, 
    const dbscan::dbPoints& hits, const larcv::Image2D& img, const int row_target, const int plane, const int row_radius ) {

    std::vector< FlashMuonTaggerAlgo::ClusterInfo_t > output;

    int rmax_boundary = row_target+row_radius;
    int rmin_boundary = row_target-row_radius;
    if ( rmax_boundary>=img.meta().rows() ) rmax_boundary = img.meta().rows()-1;
    if ( rmin_boundary<0 ) rmin_boundary = 0;

    // for each cluster we grab some info
    for ( int ic=0; ic<(int)dbscan_output.clusters.size(); ic++ ) {
      const dbscan::dbCluster&  cluster = dbscan_output.clusters.at(ic);
      // make a ClusterInfo object for this
      ClusterInfo_t info;
      info.plane = plane;
      info.cluster_idx = ic;
      info.hits_rmin_boundary = 0;
      info.hits_rmax_boundary = 0;
      info.ntime_ends = 0;
      info.row = -1;
      info.col = -1;
      info.q = -1;
      info.row_min = -1;
      info.col_min = -1;
      info.q_min = -1;
      info.row_max = -1;
      info.col_max = -1;
      info.q_max = -1;
      info.indead_region = false;
      info.npixels = (int)cluster.size();
      int max_hitidx = -1;
      for ( int ih=0; ih<cluster.size(); ih++ ) {
        int hitidx = cluster.at(ih);
        const std::vector<double>& pixel = hits.at(hitidx);
        int r = (int)pixel[1];
        int c = (int)pixel[0];
        float q = img.pixel(r,c);        
        if ( r==rmax_boundary ) {
          info.hits_rmax_boundary = 1;
          if ( info.q_max<0 || info.q_max<q ) {
            info.q_max = q;
            info.row_max = r;
            info.col_max = c;
          }
        }
        if ( r==rmin_boundary ) {
          info.hits_rmin_boundary = 1;
          if ( info.q_min<0 || info.q_min>q ) {
            info.q_min = q;
            info.row_min = r;
            info.col_min = c;
          }          
        }
	if ( r==row_target ) {
          if ( info.q<0 || q>info.q ) {
            info.q = q;
            max_hitidx = hitidx;
            info.row = r;
            info.col = c;
          }
        }
      }
      if ( info.col>=0 )
      	output.emplace_back( std::move(info) );
    }

    return output;
  }

  void FlashMuonTaggerAlgo::findPlaneMatchedClusters( const std::vector< std::vector< FlashMuonTaggerAlgo::ClusterInfo_t > >& cluster_info, 
    const std::vector<larcv::Image2D>& img_v, const float max_triarea, const std::vector<float>& z_range, const std::vector<float>& y_range, 
    std::vector< std::vector<BoundaryEndPt> >& endpts_v, std::vector< std::vector< FlashMuonTaggerAlgo::ClusterInfo_t > >& accepted_clusters ) {
    // we match clusters across plane using some markers
    //  1) consistent 2D position for the center
    //  2) if one of the t-boundaries are 2D consistent
    //  3) for cluster with no t-boundary crossings but matches in center and has a number of pixels in Y (horizontal muon) 

    // generate 3 plane combinations
    if ( fConfig.verbosity<1){
      std::cout << "  generate 3-plane intersections within "
                << " z-range=[" << z_range[0] << "," << z_range[1] << "]" 
                << " y-range=[" << y_range[0] << "," << y_range[1] << "]" 
                << std::endl;           
    }
    const int nplanes = (int)img_v.size();

    // we try all combinations of 3 planes.  factorial expansion...
    std::vector<int> numpts(nplanes);
    for (int p=0; p<nplanes; p++) {
      numpts[p] = (int)cluster_info.at(p).size();
    }
    // dont feel like writing something that generates combinatorics, going to assume 3 planes
    std::vector< std::vector<int> > combos;
    for (int p1=0; p1<numpts[0]; p1++) {
      for (int p2=0; p2<numpts[1]; p2++) {
        for ( int p3=0; p3<numpts[2]; p3++ ) {
          std::vector<int> combo(3);
          combo[0] = p1;
          combo[1] = p2;
          combo[2] = p3;
          combos.emplace_back( std::move(combo) );
        }
      }
    }

    // ok, now we go through and test all these markers
    for ( auto const& combo : combos ) {
      // test intersection of center point
      std::vector<int> wid(nplanes);
      for ( int p=0; p<nplanes; p++) {
        const ClusterInfo_t& info = cluster_info.at(p).at( combo[p] );
        wid[p] = img_v.at(p).meta().pos_x( info.col );
      }
      int crosses = 0;
      std::vector<float> poszy;
      double triangle_area;
      larcv::UBWireTool::wireIntersection( wid, poszy, triangle_area, crosses );
      if ( triangle_area > max_triarea ) continue;

      // next check the max and min boundary checks
      bool has_max_match = true;
      std::vector<int> wid_max(nplanes);
      for ( int p=0; p<nplanes; p++ ) {
        const ClusterInfo_t& info = cluster_info.at(p).at( combo[p] );
        if ( info.col_max<0 ) {
          has_max_match = false;
          break;
        }
        wid_max[p] = img_v.at(p).meta().pos_x( info.col_max );
      }
      crosses = 0;
      std::vector<float> poszy_max;
      double triarea_max;
      if ( has_max_match ) {
        larcv::UBWireTool::wireIntersection( wid_max, poszy_max, triarea_max, crosses );
        std::cout << "  test max boundary: triarea=" << triarea_max << std::endl;
        if ( crosses==0 || triarea_max>max_triarea )
          has_max_match = false;
      }

      bool has_min_match = true;
      std::vector<int> wid_min(nplanes);
      for ( int p=0; p<nplanes; p++ ) {
        const ClusterInfo_t& info = cluster_info.at(p).at( combo[p] );
        if ( info.col_min<0 ) {
          has_min_match = false;
          break;
        }
        wid_min[p] = img_v.at(p).meta().pos_x( info.col_min );
      }
      crosses = 0;
      std::vector<float> poszy_min;
      double triarea_min;
      if ( has_min_match ) {
        larcv::UBWireTool::wireIntersection( wid_min, poszy_min, triarea_min, crosses );
        std::cout << "  test min boundary: triarea=" << triarea_min << std::endl;        
        if ( crosses==0 || triarea_min>max_triarea )
          has_min_match = false;
      }

      // ok, how did we do?
      if ( !has_max_match && !has_min_match ) {
      	continue;
      }
      int nboundary_matches = 0;
      if ( has_max_match ) nboundary_matches++;
      if ( has_min_match ) nboundary_matches++;

      // we accept this combination
      if ( fConfig.verbosity<1 )
      	std::cout << "  accepted endpt: comboidx=(" << combo[0] << "," << combo[1] << "," << combo[2] << ") wires=(";
      std::vector< BoundaryEndPt > endpt_v;
      for ( int p=0; p<nplanes; p++ ) {
	const ClusterInfo_t& info = cluster_info.at(p).at( combo[p] );
        if ( fConfig.verbosity<1)
	  std::cout << img_v.at(p).meta().pos_x( info.col ) << ",";
      	BoundaryEndPt endpt( info.row, info.col );
        endpt_v.emplace_back( std::move(endpt) );
      }
      if ( fConfig.verbosity<1 )
      	std::cout << ") nboundary_matches=" << nboundary_matches << std::endl;

      endpts_v.emplace_back( std::move(endpt_v) );
    }
  }

  bool FlashMuonTaggerAlgo::findImageTrackEnds( const std::vector<larcv::Image2D>& tpc_imgs, const std::vector<larcv::Image2D>& badch_imgs,
                                                std::vector< std::vector< BoundaryEndPt > >& trackendpts ) {
    
    if ( fSearchMode!=kOutOfImage ) {
      std::cout << "[ERROR] Invalid search mode for these type of track endpoints" << std::endl;
      return false;
    }
    
    // we build fake flashes to pass into findTrackEnds
    larlite::event_opflash* faux_flashes = new larlite::event_opflash;
    std::vector<double> dummy(32,0.);
    larlite::opflash img_begin( 2400.0+36, 0, 0, 0, dummy ); // make this config pars
    larlite::opflash img_end( 2400.0+6048-36, 0, 0, 0, dummy ); // make this config pars
    faux_flashes->emplace_back( img_begin );
    faux_flashes->emplace_back( img_end );
    std::vector< larlite::event_opflash* > faux_flashes_v;
    faux_flashes_v.push_back( faux_flashes );
    bool results = flashMatchTrackEnds( faux_flashes_v, tpc_imgs, badch_imgs, trackendpts );
    delete faux_flashes;
    return results;
  }  
  

}

