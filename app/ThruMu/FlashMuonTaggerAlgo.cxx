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
	std::string modename;
	if ( fSearchMode==kAnode ) {
	  flash_tick = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
	  tick_target = flash_tick;
	  modename = "anode";
	}
	else if ( fSearchMode==kCathode ) {
	  flash_tick = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
	  tick_target = flash_tick + fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick+240.0;
	  modename = "cathode";
	}
	else if ( fSearchMode==kOutOfImage ) {
	  flash_tick  = opflash.Time();
	  tick_target = flash_tick; // dummy opflash gives first or last tick of image
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

	std::vector<float> z_range = { z_weighted+min_dist_z, z_weighted+max_dist_z }; // be more intelligent later
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
	dbscan::dbPoints* hits[nplanes];
	dbscan::dbscanOutput* cluster_info[nplanes];
	for (int p=0; p<nplanes; p++) {
	  const larcv::Image2D& img = tpc_imgs.at(p);
	  hits[p] = new dbscan::dbPoints;
	  for (int drow=-10; drow<=10; drow++) {
	    int r_current = row_target + drow;
	    if ( r_current<0 || r_current>=(int)meta.rows() )continue;
	    for (int c=0; c<(int)meta.cols(); c++) {
	      if ( img.pixel( r_current, c )>fConfig.pixel_value_threshold.at(p) ) {
		std::vector< double > pt(2);
		pt.at(0) = c;
		pt.at(1) = r_current;
		hits[p]->emplace_back(pt);
	      }
	    }
	  }
	}
	
	// we need to find clusters consistent with flashes
	
	// we find 'end points'
	std::vector< BoundaryEndPt >* endpts[3];
	for (int p=0; p<nplanes; p++) {
	  dbscan::DBSCANAlgo dbalgo;
	  cluster_info[p] = new dbscan::dbscanOutput;
	  endpts[p] = new std::vector< BoundaryEndPt >;
	  endpts[p]->clear();
	  (*cluster_info[p]) = dbalgo.scan( *(hits[p]), 5, 5.0, false, 0.0 );
	  for (int ic=0; ic<(int)cluster_info[p]->clusters.size(); ic++) {
	    BoundaryEndPt endpt;
	    bool foundend = findClusterEnds( *(cluster_info[p]), *(hits[p]), ic, row_target, p, meta, endpt );
	    if ( foundend ) endpts[p]->emplace_back( endpt );
	  }
          if (fConfig.verbosity<1 )
	    std::cout << "  flashmuontaggeralgo: plane " << p << " clusters: " << cluster_info[p]->clusters.size() << " endpoints=" << endpts[p]->size() << std::endl;
	}

	// now that we have end points for each cluster, we ID position of charge deposition in (Y,Z) and check if consistent with the flash
	std::vector< std::vector<int> > wirelists(3);
	std::vector< std::map<int,int> > wirelist_endpt_idx(3); // map from wid to endpt index
	std::vector< std::vector<float> > valid_range(2);
	valid_range[0] = z_range;
	valid_range[1] = y_range;
	for (int ip=0; ip<3; ip++) {
	  for (int ipt=0; ipt<(int)endpts[ip]->size(); ipt++) {
	    const BoundaryEndPt& endpt = endpts[ip]->at(ipt);
	    int wid = meta.pixel_width()*endpt.w;
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
		    << "score=" << areas3plane.at(ii) << std::endl;
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
// 	  for (int iw=0; iw<wire_combo.size(); iw++) {
// 	    int wid = wire_combo.at(iw);
// 	    if ( used_wid[iw].find( wid )!=used_wid[iw].end() ) {
// 	      // wire already used
// 	      valid_combo = false;
// 	      break;
// 	    }
// 	  }
	  if ( !valid_combo || combo.score>20.0 )
	    continue; // do not accept
	  
	  if ( ncombos>=max_flash_combos && combo.score>10 )
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
		int idx_endpt = endpts[missing_plane]->size();
		endpts[missing_plane]->emplace_back( std::move(endpt) );
		wirelist_endpt_idx.at(missing_plane).insert( std::pair<int,int>(wid,idx_endpt) );
	      }
	    }
	  }
	}
	

	// =================================================================================================
	// FINAL 3 PLANE COMBOS. FILL TRACK END POINTS CONTAINER WITH RESULTS
	//std::cout << "Final 3-plane combos" << std::endl;
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
	  // std::cout << "   (" << intersections3plane.at(idx).at(0) << ","
	  // 	    << intersections3plane.at(idx).at(1) << ","
	  // 	    << intersections3plane.at(idx).at(2) << ") "
	  // 	    << " endpt=(" << endptidx[0] << "," << endptidx[1] << "," << endptidx[2] << ") "
	  // 	    << " vertex (z,y)=(" << vertex3plane.at(idx).at(0) << ", " << vertex3plane.at(idx).at(1) << ")"
	  // 	    << "score=" << areas3plane.at(idx) << std::endl;
	  std::vector< BoundaryEndPt > endpt_v;
	  for (int ip=0; ip<nplanes; ip++) {
	    int eidx = endptidx.at(ip); // end point index
	    larlitecv::BoundaryEndPt endpt = endpts[ip]->at(eidx);
	    if ( fSearchMode==kAnode )
	      endpt.type = larlitecv::BoundaryEndPt::kAnode;
	    else if ( fSearchMode==kCathode )
	      endpt.type = larlitecv::BoundaryEndPt::kCathode;
	    else if ( fSearchMode==kOutOfImage ) // makes no sense here
	      endpt.type = larlitecv::BoundaryEndPt::kImageEnd;
	    endpt_v.emplace_back( std::move( endpt ) );
	  }//end of loop over planes
	  trackendpts.emplace_back( std::move( endpt_v ) );
	}
	for (int p=0; p<3; p++) {
	  delete cluster_info[p];
	  delete endpts[p];
	  delete hits[p];
	  cluster_info[p] = NULL;
	  endpts[p] = NULL;
	  hits[p] = NULL;
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
    // if ( fConfig.verbosity<2 ) {
    //   std::cout << "end points: max (r,c)=(" << tmax << ", " << wmax << ")"
    // 		<< " tw=(" << meta.pos_y(tmax) << "," << meta.pos_x(wmax) << ")"
    // 		<< "  min (r,c)=(" << tmin << "," << wmin << ")" 
    // 		<< " tw=(" << meta.pos_y(tmin) << "," << meta.pos_x(wmin) << ")"
    // 		<< "  versus: query_row=" << row_target
    // 		<< std::endl;
    // }

    // is this a marked flash-end?
    // if extrema matches the annode flash hypothesis row. mark as interesting (Score 200)
    bool success = false;
    bool tmin_isend = false;
    bool tmax_isend = false;
    if ( abs(tmin-row_target)<=fConfig.endpoint_time_neighborhood.at(plane) ) {
      //if ( fConfig.verbosity<1 ) std::cout << "MATCHED MIN END" << std::endl;
      tmin_isend = true;
    }
    if ( abs(tmax-row_target)<=fConfig.endpoint_time_neighborhood.at(plane) ) { 
      //if ( fConfig.verbosity<1 ) std::cout << "MATCHED MAX_END" << std::endl;
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

  bool FlashMuonTaggerAlgo::findImageTrackEnds( const std::vector<larcv::Image2D>& tpc_imgs, const std::vector<larcv::Image2D>& badch_imgs,
						std::vector< std::vector< BoundaryEndPt > >& trackendpts ) {
    
    if ( fSearchMode!=kOutOfImage ) {
      std::cout << "[ERROR] Invalid search mode for these type of track endpoints" << std::endl;
      return false;
    }
    
    // we build fake flashes to pass into findTrackEnds
    larlite::event_opflash* faux_flashes = new larlite::event_opflash;
    std::vector<double> dummy(32,0.);
    larlite::opflash img_begin( 2400.0+1, 0, 0, 0, dummy ); // make this config pars
    larlite::opflash img_end( 2400.0+6048-1, 0, 0, 0, dummy ); // make this config pars
    faux_flashes->emplace_back( img_begin );
    faux_flashes->emplace_back( img_end );
    std::vector< larlite::event_opflash* > faux_flashes_v;
    faux_flashes_v.push_back( faux_flashes );
    bool results = flashMatchTrackEnds( faux_flashes_v, tpc_imgs, badch_imgs, trackendpts );
    delete faux_flashes;
    return results;
  }  
  

}

