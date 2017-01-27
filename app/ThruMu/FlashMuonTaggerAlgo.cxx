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
                                                 std::vector< BoundarySpacePoint >& trackendpts ) {
    // output
    // ------
    //  trackendpts: list of vector of end pts. each (inner vector) is point on each plane.
    const int nplanes = tpc_imgs.size();
    
    if ( nplanes==0 )
      return false;
    
    if ( fConfig.verbosity<2 )
      std::cout << "Begin FlashMuonTaggerAlgo::flashMatchTrackEnds [verbsity=" << fConfig.verbosity << "]" << std::endl;

    // get a meta
    const larcv::ImageMeta& meta = tpc_imgs.at(0).meta();
    
    // loop over all the flash containers
    for ( auto& ptr_event_flash : opflashsets ) {
      // loop over flashes
      for ( auto& opflash : *ptr_event_flash ) {
        
        // determine time of flash
        float tick_target = 0;
        float flash_tick = 0;
        larlitecv::BoundaryEnd_t point_type;
        std::string modename;
        if ( fSearchMode==kAnode ) {
          flash_tick = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick + 15;
          tick_target = flash_tick;
          point_type = larlitecv::kAnode;
          modename = "anode";
        }
        else if ( fSearchMode==kCathode ) {
          flash_tick = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
          tick_target = flash_tick + fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick+240.0;
          point_type = larlitecv::kCathode;          
          modename = "cathode";
        }
        else if ( fSearchMode==kOutOfImage ) {
          flash_tick  = opflash.Time();
          tick_target = flash_tick; // dummy opflash gives first or last tick of image
          point_type = larlitecv::kImageEnd;          
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
        
        if ( fConfig.verbosity<2 ) {
          std::cout << "============================================================================================" << std::endl;
          std::cout << "[Opflash search] " << modename << " mode" << std::endl;
          std::cout << "  opflash: "
                    << " flash_tick=" << flash_tick
                    << " tick_target=" << tick_target
                    << " row_target=" << row_target
                    << " qtot= " << qtot
                    << " z_range=[" << z_range[0] << "," << z_range[1] << "] "
                    << " w_range=[" << int(z_range[0]/0.3/meta.pixel_width()) << "," << int(z_range[1]/0.3/meta.pixel_width()) << "] "
                    << " drift_t=" << fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick+200.0 << " ticks"
                    << std::endl;
        }
        
        // first we make clusters in the neighborhood of the target_row
        std::vector<dbscan::dbPoints> hits;
        for (int p=0; p<nplanes; p++) {
          const larcv::Image2D& img = tpc_imgs.at(p);
          dbscan::dbPoints planehits;
          for (int drow=-20; drow<=20; drow++) {
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
        
        // on each plane, separately, we cluster hits consistent with the flash
        std::vector<dbscan::dbscanOutput> dbscan_output; 
        std::vector< std::vector<ClusterInfo_t> > cluster_info;
        for (int p=0; p<nplanes; p++) {
          dbscan::DBSCANAlgo dbalgo;
          dbscan::dbscanOutput plane_dbscan_output      = dbalgo.scan( hits.at(p), 5, 5.0, false, 0.0 );
          std::vector<ClusterInfo_t> plane_cluster_info = analyzeClusters( plane_dbscan_output, hits.at(p), tpc_imgs.at(p), row_target, p, 5 );
          for ( auto& info : plane_cluster_info ) {
            if ( fConfig.verbosity<1) {
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
          }
          dbscan_output.emplace_back( std::move(plane_dbscan_output) );
          cluster_info.emplace_back( std::move(plane_cluster_info) );
        }

        // now we match up clusters from across the 3 planes
        std::vector< BoundarySpacePoint > endpts_v;
        std::vector< std::vector<ClusterInfo_t> > accepted_cluster_matches;        
        findPlaneMatchedClusters( cluster_info, tpc_imgs, badch_imgs, fConfig.max_triarea, z_range, y_range, endpts_v, accepted_cluster_matches );
        if ( fConfig.verbosity<2 ) {
          std::cout << "number of plane-matched clusters: " << accepted_cluster_matches.size() << std::endl;
          std::cout << "number of endpts: " << endpts_v.size() << std::endl;
        }

        // now we try to remove false positives
        std::vector< int > passing_clusters;                
        filterClusters( accepted_cluster_matches, tpc_imgs, 20, 20, 20, passing_clusters );

        // transfer to output container
        if ( fConfig.verbosity<2 )
          std::cout << "Generated End Points" << std::endl;
        for (int idx_cluster=0; idx_cluster<(int)passing_clusters.size(); idx_cluster++ ) {
          if ( passing_clusters.at(idx_cluster)==0 ) continue;
          BoundarySpacePoint& intersection = endpts_v.at(idx_cluster);

          std::cout << "  (" << tpc_imgs.at(0).meta().pos_x(intersection[0].col) 
            << "," << tpc_imgs.at(1).meta().pos_x(intersection[1].col)
            << "," << tpc_imgs.at(2).meta().pos_x(intersection[2].col) << ")" 
            << std::endl;
          for ( auto& plane_pt : intersection )
            plane_pt.type = point_type;
          trackendpts.emplace_back( std::move( intersection ) );
        }

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
  
  BoundaryEnd_t FlashMuonTaggerAlgo::SearchModeToEndType( SearchMode_t mode ) {
    switch(mode) {
      case FlashMuonTaggerAlgo::kAnode:
        return larlitecv::kAnode;
        break;
      case FlashMuonTaggerAlgo::kCathode:
        return larlitecv::kCathode;
        break;
      case FlashMuonTaggerAlgo::kOutOfImage:
        return larlitecv::kImageEnd;
        break;
      default:
        throw std::runtime_error("FlashMuonTaggerAlgo::SearchModeToEndType unknown search mode.");
        break;
    }
    return larlitecv::kUndefined;
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
      endpt.row = tmin;
      endpt.col = wmin;
      endpt.type = SearchModeToEndType( fSearchMode );
    }
    else if ( (tmax_isend && !tmin_isend) || (tmax_isend && tmin_isend && abs(tmax-row_target)<abs(tmin-row_target)) ) {
      success = true;
      endpt.row = tmax;
      endpt.col = wmax;
      endpt.type = SearchModeToEndType( fSearchMode );
    }      
        
    return success;
  }

  void FlashMuonTaggerAlgo::generate3PlaneIntersections(  const std::vector< std::vector< ClusterInfo_t > >& cluster_info, const std::vector<larcv::Image2D>& img_v, 
    const std::vector<float>& z_range, const std::vector<float>& y_range, const float max_triarea,
    std::vector< BoundarySpacePoint >& endpts_v, std::vector< std::vector< ClusterInfo_t > >& accepted_clusters,
    std::vector< std::vector<int> >& cluster_used ) {
    // input
    //  vector< vector<cluster_info> >:  list of cluster info for each plane
    //  vector< image2d >: list of images, 1 for each plane
    // output
    //  vector< vector<BoundaryEndPt> > : list of plane-matched end-pts. inner vector has 3 endpts each, one for each plane.
    //  vector< vector<ClusterInfo_t> > : list of plane-matched cluster info. mirrors above. inner vector has 3 info objects each, one for each plane

     // generate 3 plane combinations
    if ( fConfig.verbosity<1){
      std::cout << " ** generate 3-plane intersections within "
                << " z-range=[" << z_range[0] << "," << z_range[1] << "]" 
                << " y-range=[" << y_range[0] << "," << y_range[1] << "] **" 
                << std::endl;           
    }
    const int nplanes = (int)img_v.size();

    // we try all combinations of 3 planes.  factorial expansion...
    std::vector<int> numpts(nplanes);
    cluster_used.clear();
    cluster_used.resize(nplanes);
    for (int p=0; p<nplanes; p++) {
      numpts[p] = (int)cluster_info.at(p).size();
      cluster_used[p].resize( numpts[p], 0);
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
      //if ( has_max_match ) {
      //  larcv::UBWireTool::wireIntersection( wid_max, poszy_max, triarea_max, crosses );
      //  if ( crosses==0 || triarea_max>max_triarea )
      //    has_max_match = false;
      //}

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
      //if ( has_min_match ) {
      //  larcv::UBWireTool::wireIntersection( wid_min, poszy_min, triarea_min, crosses );
      //  if ( crosses==0 || triarea_min>max_triarea )
      //    has_min_match = false;
      //}

      // ok, how did we do?
      if ( !has_max_match && !has_min_match ) {
        continue;
      }
      int nboundary_matches = 0;
      if ( has_max_match ) nboundary_matches++;
      if ( has_min_match ) nboundary_matches++;

      // we accept this combination
      if ( fConfig.verbosity<1 )
        std::cout << "  accepted 3-plane endpt: comboidx=(" << combo[0] << "," << combo[1] << "," << combo[2] << ") wires=(";
      std::vector< BoundaryEndPt > endpt_v;
      std::vector< ClusterInfo_t > info_copy_v;
      for ( int p=0; p<nplanes; p++ ) {
        const ClusterInfo_t& info = cluster_info.at(p).at( combo[p] );
        cluster_used[p][combo[p]] = 1;
        if ( fConfig.verbosity<1)
          std::cout << img_v.at(p).meta().pos_x( info.col ) << ",";
        BoundaryEndPt endpt( info.row, info.col, SearchModeToEndType(fSearchMode) );
        endpt_v.emplace_back( std::move(endpt) );
        info_copy_v.push_back( info );
      }

      BoundarySpacePoint spacepoint( SearchModeToEndType(fSearchMode), std::move(endpt_v) );


      if ( fConfig.verbosity<1 )
        std::cout << ") nboundary_matches=" << nboundary_matches << " max=" << triarea_max <<  " min=" << triarea_min << std::endl;

      endpts_v.emplace_back( std::move(spacepoint) );
      accepted_clusters.emplace_back( std::move(info_copy_v) );
    }
  }

  
  void FlashMuonTaggerAlgo::generate2PlaneIntersections( const std::vector< std::vector< ClusterInfo_t > >& cluster_info, 
    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
    const std::vector<float>& z_range, const std::vector<float>& y_range,
    std::vector< BoundarySpacePoint >& endpts_v, std::vector< std::vector< ClusterInfo_t > >& accepted_clusters,
    std::vector< std::vector<int> >& cluster_used ) {

    if ( fConfig.verbosity<1){
      std::cout << "  generate 3-plane intersections within "
                << " z-range=[" << z_range[0] << "," << z_range[1] << "]" 
                << " y-range=[" << y_range[0] << "," << y_range[1] << "]" 
                << std::endl;           
    }

    // we try all combinations of 2 planes.  factorial expansion...
    const int nplanes = img_v.size();
    std::vector<int> numpts(nplanes);
    for (int p=0; p<nplanes; p++) {
      numpts[p] = (int)cluster_info.at(p).size();
    }
    // find 2 plane combinations + dead ch
    std::vector< std::vector<int> > combos;

    std::vector< std::vector<ClusterInfo_t> > proposed_clusters;
    std::vector< BoundarySpacePoint > proposed_endpts;

    for (int p1=0; p1<nplanes; p1++) {
      for (int p2=p1+1; p2<nplanes; p2++) {

        for ( int endpt1=0; endpt1<numpts[p1]; endpt1++) {
          // check if this endpt is already in use
          //if ( cluster_used[p1][endpt1]==1 ) continue;

          for ( int endpt2=0; endpt2<numpts[p2]; endpt2++ ) {
            // check if this is also in use
            //if ( cluster_used[p2][endpt2]==1 ) continue;

            // OK, lets see what wire we have
            int wid1 = img_v.at(p1).meta().pos_x( cluster_info[p1][endpt1].col );
            int wid2 = img_v.at(p2).meta().pos_x( cluster_info[p2][endpt2].col );

            std::vector<float> poszy;
            int crosses;
            int otherplane;
            int otherwire;
            larcv::UBWireTool::getMissingWireAndPlane( p1, wid1, p2, wid2, otherplane, otherwire, poszy, crosses );

            int badch_state = 0;
            if ( crosses>0 && otherwire>=0 && otherplane>=0 && otherwire<img_v.at(otherplane).meta().max_x() )
              badch_state = badch_v.at(otherplane).pixel( cluster_info[p1][endpt1].row, img_v.at(otherplane).meta().col(otherwire) );

            if (fConfig.verbosity<1 ) {
              std::cout << "p1=" << p1 << " wid1=" << wid1
                << " p2=" << p2 << " wid2=" << wid2 
                << " other plane=" << otherplane << " other wire=" << otherwire << " badch=" << badch_state
                << std::endl;
            }

            // is there an intersection?
            if ( crosses==0 )
              continue;

            // is this wire in a badch region?
            if ( badch_state== 0 ) {
              // not a badch move on
              continue;
            }

            // we test the good planes to see if they have also have a match away from the flash point


            // if so, define the combo
            ClusterInfo_t badch_info;
            badch_info.plane = otherplane;
            badch_info.indead_region = true;
            badch_info.row = cluster_info[p1][endpt1].row;
            badch_info.col = img_v.at(otherplane).meta().col(otherwire);
            std::vector< ClusterInfo_t > info_v(nplanes);
            info_v[p1] = cluster_info[p1][endpt1];
            info_v[p2] = cluster_info[p2][endpt2];
            info_v[otherplane] = badch_info;
            cluster_used[p1][endpt1] = 1;
            cluster_used[p2][endpt2] = 1;

            std::vector< BoundaryEndPt > endpts;
            for ( int p=0; p<nplanes; p++ ) {
              const ClusterInfo_t& info = info_v.at(p);
              BoundaryEndPt endpt( info.row, info.col, SearchModeToEndType(fSearchMode) );
              endpts.emplace_back( std::move(endpt) );
            }

            BoundarySpacePoint spacepoint( SearchModeToEndType(fSearchMode), std::move(endpts) );

            proposed_endpts.emplace_back( std::move(spacepoint) );
            proposed_clusters.emplace_back( std::move(info_v) );


            break; // we already matched endpt1, so break this loop
          }//end of loop over endpt2
        }
      }
    }

    for ( int idx=0; idx<(int)proposed_clusters.size(); idx++ ) {
      // for the clusters that was found, we do the same type of check as the 3-plane matches
      accepted_clusters.emplace_back( std::move(proposed_clusters.at(idx)) );
      endpts_v.emplace_back( std::move(proposed_endpts.at(idx)) );
    }

  }

  std::vector<FlashMuonTaggerAlgo::ClusterInfo_t> FlashMuonTaggerAlgo::analyzeClusters( const dbscan::dbscanOutput& dbscan_output, 
    const dbscan::dbPoints& hits, const larcv::Image2D& img, const int row_target, const int plane, const int row_radius ) {

    std::vector< FlashMuonTaggerAlgo::ClusterInfo_t > output;

    int rmax_boundary = row_target+row_radius;
    int rmin_boundary = row_target-row_radius;
    if ( rmax_boundary>=(int)img.meta().rows() ) rmax_boundary = (int)img.meta().rows()-1;
    if ( rmin_boundary<0 ) rmin_boundary = 0;

    // for each cluster we grab some info
    //std::cout << "number of cluster on plane=" << plane << ": " << dbscan_output.clusters.size() << std::endl;
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
      for ( int ih=0; ih<(int)cluster.size(); ih++ ) {
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
          if ( info.q_min<0 || info.q_min<q ) {
            info.q_min = q;
            info.row_min = r;
            info.col_min = c;
          }          
        }
        if ( r==row_target ) {
          if ( info.q<0 || info.q<q ) {
            info.q = q;
            max_hitidx = hitidx;
            info.row = r;
            info.col = c;
          }
        }
      }
      if ( max_hitidx>=0 )
        output.emplace_back( std::move(info) );
    }

    return output;
  }

  void FlashMuonTaggerAlgo::findPlaneMatchedClusters( const std::vector< std::vector< FlashMuonTaggerAlgo::ClusterInfo_t > >& cluster_info, 
    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
    const float max_triarea, const std::vector<float>& z_range, const std::vector<float>& y_range, 
    std::vector< BoundarySpacePoint >& endpts_v, std::vector< std::vector< FlashMuonTaggerAlgo::ClusterInfo_t > >& accepted_clusters ) {
    // we match clusters across plane using some markers
    //  1) consistent 2D position for the center
    //  2) if one of the t-boundaries are 2D consistent
    //  3) for cluster with no t-boundary crossings but matches in center and has a number of pixels in Y (horizontal muon) 

    // generate 3 plane combinations
    std::vector< std::vector<int> > cluster_used_for3planematch;
    generate3PlaneIntersections( cluster_info, img_v, z_range, y_range, max_triarea, endpts_v, accepted_clusters, cluster_used_for3planematch );

    // if remaining clusters, look for 2 plane combinations + dead channels
    generate2PlaneIntersections( cluster_info, img_v, badch_v, z_range, y_range, endpts_v, accepted_clusters, cluster_used_for3planematch );
  }

  void FlashMuonTaggerAlgo::filterClusters( const std::vector< std::vector< ClusterInfo_t > >& accepted_cluster_matches, 
    const std::vector<larcv::Image2D>& img_v, const int rmax_window, const int rmin_window, const int col_width,
    std::vector< int >& cluster_passed ) {
    // we want to remove false positive flash-end detections.
    // this means removing flash-tags that occur in the middle of a track.
    // such tags have the consequence of causing many multiple tracks, greatly inflatingn the number of aStar searchs
    // (though in principle, these repeate 3D tracks can be removed/merged -- but let's try to get rid of them here)
    cluster_passed.resize( accepted_cluster_matches.size(), 0 );

    int idx_cluster = 0;
    for ( auto const& info_v : accepted_cluster_matches ) {
      // to determine if track end, we do something very simple:
      // we define a window around the cluster pt in all planes
      // within the window we collect charge and cluster.
      // we find the extrema of all the points
      // the extrema tells us whether to check the horizontally of vertically
      // does both or only one extrema reach the end?
      // if only one, it passes
      // if both, it fails as being midtrack

      const int nplanes = (int)img_v.size();
      int row = info_v.front().row; // clusters on each plane should have same row
      // define the window
      int row_start = row - rmin_window;
      int row_end   = row + rmax_window;
      if ( row_start<0 ) row_start = 0;
      if ( row_end>=(int)img_v.at(0).meta().rows() ) row_end = (int)img_v.at(0).meta().rows()-1;

      dbscan::DBSCANAlgo algo;
      
      // collect hits in neighborhood
      std::vector<int> boundaries_reached(nplanes,0);
      int reached_tick_top = 0;
      int reached_tick_bot = 0;

      for (int p=0; p<nplanes;p++) {

        if ( info_v.at(p).indead_region==true )
          continue;

        dbscan::dbPoints planepts;
        int col = info_v.at(p).col;
        int col_start = col-col_width;
        int col_end   = col+col_width;
        if ( col_start<0 ) col_start=0;
        if ( col_end>=(int)img_v.at(p).meta().cols()) col_end = (int)img_v.at(p).meta().cols()-1;

        for ( int r=row_start; r<=row_end; r++ ) {
          if ( r<0 || r>=(int)img_v.at(p).meta().rows() ) continue;
          for (int c=col_start; c<=col_end; c++ ) {
            if ( c<0 || c>=(int)img_v.at(p).meta().cols() ) continue;

            if ( img_v.at(p).pixel(r,c)>fConfig.pixel_value_threshold.at(p) ) {
              std::vector<double> hit(2);
              hit[0] = c;
              hit[1] = r;
              planepts.emplace_back( std::move(hit) );
            }

          }
        }
        dbscan::dbscanOutput clout = algo.scan( planepts, 3, 5.0, false, 0.0 );
        std::vector<double> centerpt(2);
        centerpt[0] = col;
        centerpt[1] = row;
        int matching_cluster = clout.findMatchingCluster( centerpt, planepts, 3.0 );
        if ( matching_cluster<0 )
          continue;

        int cluster_size = clout.clusters.at(matching_cluster).size();
        if ( cluster_size==0 ) {
          throw std::runtime_error("previous established cluster now missing?");
        }

        // initialize extrema search
        dbscan::ClusterExtrema extrema = dbscan::ClusterExtrema::FindClusterExtrema( clout.clusters.at(matching_cluster), planepts );


        int num_boundaries_reached = 0;
        if ( extrema.bottommost()[1]<=row_start ) {
          num_boundaries_reached++;
          reached_tick_top++;
        }
        else if ( (extrema.leftmost()[0]<=col_start && extrema.leftmost()[1]<row && extrema.leftmost()[1]>row_start ) 
          || (extrema.rightmost()[0]>=col_end && extrema.rightmost()[1]<row && extrema.rightmost()[1]>row_start) ) {
          num_boundaries_reached++;
          reached_tick_top++;
        }

        if ( extrema.topmost()[1]>=row_end   ) {
          num_boundaries_reached++;
          reached_tick_bot++;
        }
        else if ( (extrema.leftmost()[0]<=col_start && extrema.leftmost()[1]>row && extrema.rightmost()[1]<row_end) 
          || (extrema.rightmost()[0]>=col_end && extrema.rightmost()[1]>row && extrema.rightmost()[1]<row_end) ) {
          num_boundaries_reached++;
          reached_tick_bot++;
        }

        boundaries_reached[p] = num_boundaries_reached;
        if ( fConfig.verbosity<1 ) {
          std::cout << "  plane=" << p << " tick-top extreme=" << img_v.at(p).meta().pos_y(extrema.bottommost()[1]) << " <= row_top=" << img_v.at(p).meta().pos_y(row_start) << std::endl;
          std::cout << "  plane=" << p << " tick-bot extreme=" << img_v.at(p).meta().pos_y(extrema.topmost()[1]) << " >= row_bot=" << img_v.at(p).meta().pos_y(row_end) << std::endl;        
          std::cout << "  plane=" << p << " tick-lft extreme=" << img_v.at(p).meta().pos_y(extrema.leftmost()[1]) 
            << " ; col=" << img_v.at(p).meta().pos_x(extrema.leftmost()[0]) << " >= " << col_start << std::endl;
          std::cout << "  plane=" << p << " tick-rgt extreme=" << img_v.at(p).meta().pos_y(extrema.rightmost()[1]) 
            << " ; col=" << img_v.at(p).meta().pos_x(extrema.rightmost()[0]) << " <= " << col_end << std::endl;        
        }

      }//end of plane loop


      int nplanes_with_crossing_cluster = 0;
      for (int p=0; p<nplanes; p++) {
        if ( boundaries_reached[p]>=2 ) nplanes_with_crossing_cluster++;
      }

      if ( nplanes_with_crossing_cluster<=1 || fSearchMode==kOutOfImage ) {
        // conditions for passing
        if ( (fSearchMode==kCathode || (fSearchMode==kOutOfImage && fOutOfImageMode==kFront) ) && reached_tick_bot>=2 )
          cluster_passed[idx_cluster] = 1;
        else if ( (fSearchMode==kAnode || (fSearchMode==kOutOfImage && fOutOfImageMode==kBack)) && reached_tick_top>=2 )
          cluster_passed[idx_cluster] = 1;
      }
      if ( fConfig.verbosity<1 ) {
        std::cout << "filter cluster#" << idx_cluster << " tick=" << img_v.at(0).meta().pos_y(row) 
          << " cols=(" << info_v[0].col << "," << info_v[1].col << "," << info_v[2].col << ")"
          << " [planes with crossing cluster]=" << nplanes_with_crossing_cluster
          << " [planes reached top tick]=" << reached_tick_top
          << " [planes reached bot tick]=" << reached_tick_bot
          << " [cluster_passed]=" << cluster_passed[idx_cluster]
          << std::endl;
      }
      idx_cluster++;
    }//end of cluster loop

  }// end of function

  bool FlashMuonTaggerAlgo::findImageTrackEnds( const std::vector<larcv::Image2D>& tpc_imgs, const std::vector<larcv::Image2D>& badch_imgs,
                                                std::vector<BoundarySpacePoint >& trackendpts ) {
    
    if ( fSearchMode!=kOutOfImage ) {
      std::cout << "[ERROR] Invalid search mode for these type of track endpoints" << std::endl;
      return false;
    }
    
    // we build fake flashes to pass into findTrackEnds
    
    std::vector<double> dummy(32,0.);

    larlite::event_opflash* faux_flash_front = new larlite::event_opflash;
    larlite::opflash img_begin( 2415.0, 0, 0, 0, dummy );    
    faux_flash_front->emplace_back( img_begin );
    std::vector< larlite::event_opflash* > faux_flash_front_v;    
    faux_flash_front_v.push_back( faux_flash_front );
    fOutOfImageMode = kFront;
    bool results = flashMatchTrackEnds( faux_flash_front_v, tpc_imgs, badch_imgs, trackendpts );

    larlite::event_opflash* faux_flash_back = new larlite::event_opflash;        
    larlite::opflash img_end( 2400.0+6048-36, 0, 0, 0, dummy ); // make this config pars
    faux_flash_back->emplace_back( img_end );
    std::vector< larlite::event_opflash* > faux_flash_back_v;    
    faux_flash_back_v.push_back( faux_flash_back );
    fOutOfImageMode = kBack;
    results = flashMatchTrackEnds( faux_flash_back_v, tpc_imgs, badch_imgs, trackendpts );

    delete faux_flash_back;
    delete faux_flash_front;

    return results;
  }  
  

}

