#include "FlashMuonTaggerAlgo.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"

namespace larlitecv {
  
  bool FlashMuonTaggerAlgo::findTrackEnds( const std::vector< larlite::event_opflash* >& opflashsets, const larcv::Image2D& tpc_img,
					   std::vector< BoundaryEndPt >& trackendpts, larcv::Image2D& markedimg ) {
    
    markedimg = std::move( larcv::Image2D(tpc_img.meta()) );

    for ( auto& ptr_event_flash : opflashsets ) {
      for ( auto& opflash : *ptr_event_flash ) {
	
	const larcv::ImageMeta& meta = tpc_img.meta();
	int plane = meta.plane();
	
	float tick_target = 0;
	std::string modename;
	if ( fSearchMode==kAnode ) {
	  tick_target = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
	  modename = "anode";
	}
	else if ( fSearchMode==kCathode ) {
	  tick_target = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
	  tick_target += fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick + 240.0;
	  modename = "cathode";
	}
	else if ( fSearchMode==kOutOfImage ) {
	  tick_target = opflash.Time(); // dummy opflash gives first or last tick of image
	  modename = "image end";
	}
	else {
	  std::cout << "[ERROR] wrong search mode" << std::endl;
	  return false;
	}

	// check if we can search for this opflash
	if ( tick_target<meta.min_y() || tick_target>=meta.max_y() )
	  continue;

	int row_target = meta.row( tick_target );

	if ( fConfig.verbosity<1 ) {
	  std::cout << "============================================================================================" << std::endl;
	  std::cout << "[Opflash search] " << modename << " mode" << std::endl;
	  std::cout << "  opflash: p= " << plane 
		    << " tick_target=" << tick_target
		    << " row_target=" << row_target
		    << " drift_t=" << fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick << " ticks"
		    << std::endl;
	}
	
	// scan across wires, looking for a hit
	for (int iwire=0; iwire<meta.cols(); iwire++) {
	  if ( markedimg.pixel( row_target, iwire )<= 1.0 // not visited
	       && tpc_img.pixel( row_target, iwire )>fConfig.pixel_value_threshold[plane] ) {// above threshold  
	    
	    // cluster and return cluster that contains the query point
	    dbscan::dbscanOutput clout;
	    dbscan::dbPoints winpoints;
	    int containing_cluster = -1;
	    bool ok = getClusters( tpc_img, row_target, iwire, clout, winpoints, containing_cluster );
	    
	    if ( !ok ) continue;

	    // find extremities of cluster (in time)
	    int tmax = -1;
	    int tmin = -1;
	    int wmax = -1;
	    int wmin = -1;
	    for ( int ichit=0; ichit<clout.clusters.at(containing_cluster).size(); ichit++ ) {
	      int hitidx = clout.clusters.at(containing_cluster).at(ichit);
	      int x_ = (int)winpoints.at(hitidx).at(0)+0.1;
	      int y_ = (int)winpoints.at(hitidx).at(1)+0.1;
	      //markedimg.set_pixel( y_, x_, 10.0 );
	      if ( tmax==-1 || y_>tmax ) { tmax = y_; wmax = x_; };
	      if ( tmin==-1 || y_<tmin ) { tmin = y_; wmin = x_; };
	    }
	    if ( fConfig.verbosity<2 ) {
	      std::cout << "end points: max (r,c)=(" << tmax << ", " << wmax << ")"
			<< " tw=(" << meta.pos_y(tmax) << "," << meta.pos_x(wmax) << ")"
			<< "  min (r,c)=(" << tmin << "," << wmin << ")" 
			<< " tw=(" << meta.pos_y(tmin) << "," << meta.pos_x(wmin) << ")"
			<< "  versus: query_row=" << row_target
			<< std::endl;
	    }

	    // is this a marked flash-end?
	    // if extrema matches the annode flash hypothesis row. mark as interesting (Score 200)
	    bool success = false;
	    if ( abs(tmin-row_target)<=fConfig.endpoint_time_neighborhood.at(plane) ) {
	      if ( fConfig.verbosity<1 ) std::cout << "MATCHED MIN END" << std::endl;
	      markedimg.set_pixel( tmin, wmin, 200.0 );
	      success = true;
	      BoundaryEndPt endpt;
	      endpt.t = tmin;
	      endpt.w = wmin;
	      if ( fSearchMode==kAnode ) endpt.type = BoundaryEndPt::kAnode;
	      else if ( fSearchMode==kCathode ) endpt.type = BoundaryEndPt::kCathode;
	      else if ( fSearchMode==kOutOfImage ) endpt.type = BoundaryEndPt::kImageEnd;
	      trackendpts.emplace_back( std::move(endpt) );
	    }
	    else if ( abs(tmax-row_target)<=fConfig.endpoint_time_neighborhood.at(plane) ) { 
	      if ( fConfig.verbosity<1 ) std::cout << "MATCHED MAX_END" << std::endl;
	      success =true;
	      markedimg.set_pixel( tmax, wmax, 200.0 );
	      BoundaryEndPt endpt;
	      endpt.t = tmax;
	      endpt.w = wmax;
	      if ( fSearchMode==kAnode ) endpt.type = BoundaryEndPt::kAnode;
	      else if ( fSearchMode==kCathode ) endpt.type = BoundaryEndPt::kCathode;
	      else if ( fSearchMode==kOutOfImage ) endpt.type = BoundaryEndPt::kImageEnd;
	      trackendpts.emplace_back( std::move(endpt) );
	    }
	    
	    if ( success ) {
	      // mark good cluster, so we don't use it again
	      for ( int ichit=0; ichit<clout.clusters.at(containing_cluster).size(); ichit++ ) {
		int hitidx = clout.clusters.at(containing_cluster).at(ichit);
		int x_ = (int)winpoints.at(hitidx).at(0)+0.1;
		int y_ = (int)winpoints.at(hitidx).at(1)+0.1;
		if ( markedimg.pixel( y_, x_ )<100 )
		  markedimg.set_pixel( y_, x_, 10.0 );      
	      }	      
	    }
	    
	  }//end of if point is interesting
	}//end of loop over wires
      }//end of loop over flashes
    }//end of loop over flash containers
    return true;
  }//end of marked output
  
  
  bool FlashMuonTaggerAlgo::flashMatchTrackEnds( const std::vector< larlite::event_opflash* >& opflashsets, const std::vector<larcv::Image2D>& tpc_imgs,
						 std::vector< std::vector< BoundaryEndPt > >& trackendpts, std::vector< larcv::Image2D >& markedimgs ) {

    const int nplanes = tpc_imgs.size();
    trackendpts.clear();
    trackendpts.resize(nplanes);
    
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
	  std::cout << "[ERROR] wrong search mode" << std::endl;
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
	
	int row_target = meta.row( tick_target );
	
	if ( true || fConfig.verbosity<1 ) {
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
	    if ( r_current<0 || r_current>=meta.rows() )continue;
	    for (int c=0; c<meta.cols(); c++) {
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
	  (*cluster_info[p]) = dbalgo.scan( *(hits[p]), 5, 5.0, false, 0.0 );
	  for (int ic=0; ic<cluster_info[p]->clusters.size(); ic++) {
	    BoundaryEndPt endpt;
	    bool foundend = findClusterEnds( *(cluster_info[p]), *(hits[p]), ic, row_target, p, meta, endpt, markedimgs.at(p) );
	    if ( foundend ) endpts[p]->emplace_back( endpt );
	  }
	  std::cout << "  flashmuontaggeralgo: plane " << p << " clusters: " << cluster_info[p]->clusters.size() << " endpoints=" << endpts[0]->size() << std::endl;
	}
	
	// now that we have end points for each cluster, we ID position of charge deposition in (Y,Z) and check if consistent with the flash
	std::vector< std::vector<int> > wirelists(3);
	std::vector< std::map<int,int> > wirelist_endpt_idx(3); // map from wid to endpt index
	std::vector< std::vector<float> > valid_range(2);
	valid_range[0] = z_range;
	valid_range[1] = y_range;
	for (int ip=0; ip<3; ip++) {
	  for (int ipt=0; ipt<endpts[ip]->size(); ipt++) {
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
	
	findWireIntersections( wirelists, valid_range, intersections3plane, vertex3plane, areas3plane, intersections2plane, vertex2plane );

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
	
	// now we go through and find our best match(es)
	int max_flash_combos = 3;
	int ncombos = 0;
	std::vector< int > final_combos_idx;

	// go through 3 plane data
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
	  for (int iw=0; iw<wire_combo.size(); iw++) {
	    int wid = wire_combo.at(iw);
	    used_wid[iw].insert( wid );
	  }
	  
	  final_combos_idx.push_back( combo.idx ); // copy combo
	  ncombos++;
	}
	
	// final 3-plane combos. fill the trackendpts container
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
	  std::cout << "   (" << intersections3plane.at(idx).at(0) << ","
		    << intersections3plane.at(idx).at(1) << ","
		    << intersections3plane.at(idx).at(2) << ") "
		    << " endpt=(" << endptidx[0] << "," << endptidx[1] << "," << endptidx[2] << ") "
		    << " vertex (z,y)=(" << vertex3plane.at(idx).at(0) << ", " << vertex3plane.at(idx).at(1) << ")"
		    << "score=" << areas3plane.at(idx) << std::endl;
	  for (int ip=0; ip<nplanes; ip++) {
	    int wid = intersections3plane.at(idx).at(ip);
	    int eidx = endptidx.at(ip);
	    larlitecv::BoundaryEndPt endpt = endpts[ip]->at(eidx);
	    if ( fSearchMode==kAnode )
	      endpt.type = larlitecv::BoundaryEndPt::kAnode;
	    else if ( fSearchMode==kCathode )
	      endpt.type = larlitecv::BoundaryEndPt::kCathode;
	    else if ( fSearchMode==kOutOfImage ) // makes no sense here
	      endpt.type = larlitecv::BoundaryEndPt::kImageEnd;
	    trackendpts.at( ip ).emplace_back( endpt );
	  }
	}

      }//end of opflashes loop
    }//end of opflashsets    
    return true;
  }
  
  
  // subalgos
  bool FlashMuonTaggerAlgo::getClusters( const larcv::Image2D& tpc_img, int query_row, int query_col, 
					 dbscan::dbscanOutput& cluster_info,  dbscan::dbPoints& winpoints, int& containing_cluster) {
    winpoints.clear();

    const larcv::ImageMeta& meta = tpc_img.meta();
    int plane = meta.plane();
    
    // new, unexplored region!
    // we define a window around this point: (query_row, c)
    float t1 = meta.pos_y(query_row);
    float w  = meta.pos_x(query_col);

    if ( fConfig.verbosity<2) 
      std::cout << "[flashmuontagger: clustering region, return cluster and cluster endpoints]" << std::endl;

    int centeridx = -1;
    for (int dwire=-fConfig.clustering_wire_neighborhood.at(plane); dwire<=fConfig.clustering_wire_neighborhood.at(plane); dwire++) {
      for (int dtwin=-fConfig.clustering_time_neighborhood.at(plane); dtwin<=fConfig.clustering_time_neighborhood.at(plane); dtwin++) {
	int r = query_row+dtwin;
	int c = query_col+dwire;
	// check validity
	if ( r>=0 && r<meta.rows() && c>=0 && c<meta.cols() && tpc_img.pixel(r,c)>fConfig.pixel_value_threshold[plane] ) {
	  // mark as visited
	  //stage1_annode_hits.at(plane).set_pixel( r, c, 10.0 );
	  std::vector<double> pt(2,0.0);
	  pt[0] = c;
	  pt[1] = r;
	  winpoints.emplace_back( pt );
	  if ( dtwin==0 && dwire==0 ) centeridx = winpoints.size()-1;
	}// if valid point
      }//end of time window
    }//end of wire window
    if ( fConfig.verbosity<1 ) {
      std::cout << "  exploring around: plane=" << plane
		<< " (c,r)=" << query_col << ", " << query_row << " (w,t)=(" << w << "," << t1 << ")  npoints=" << winpoints.size() << std::endl;
    }

    // clustering
    dbscan::DBSCANAlgo dbalgo;
    cluster_info = dbalgo.scan( winpoints, 3, 3.0, false, 0.0 );

    // which cluster is the center point in?
    containing_cluster = cluster_info.clusterid.at(centeridx);
    if ( fConfig.verbosity<2 ) {
      std::cout << "  connected clusterid=" << containing_cluster  << std::endl;
      for (int ic=0; ic<cluster_info.clusters.size(); ic++) {
	std::cout << "    clusterid=" << ic << ", size=" << cluster_info.clusters.at(ic).size() << std::endl;
      }
    }
		  
    // cluster that hit is a part of is too small
    if ( containing_cluster<0 || cluster_info.clusters.at(containing_cluster).size()<5 )
      return false;

    if ( fConfig.verbosity<1 )
      std::cout << "  valid connected cluster: " << containing_cluster << " cluster size=" << cluster_info.clusters.at(containing_cluster).size() << std::endl;
    
    return true;
    
      

  }


  bool FlashMuonTaggerAlgo::findImageBoundaryEnds( const larcv::Image2D& tpc_img, std::vector< BoundaryEndPt >& trackendpts, larcv::Image2D& markedimg ) {
						   
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
    bool result = findTrackEnds( faux_flashes_v, tpc_img, trackendpts, markedimg );
    delete faux_flashes;
    return result;
  }  

  bool FlashMuonTaggerAlgo::findImageTrackEnds( const std::vector<larcv::Image2D>& tpc_imgs, std::vector< std::vector< BoundaryEndPt > >& trackendpts, std::vector< larcv::Image2D >& markedimgs ) {
						   
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
    //bool result = findTrackEnds( faux_flashes_v, tpc_img, trackendpts, markedimg );
    bool results = flashMatchTrackEnds( faux_flashes_v, tpc_imgs, trackendpts, markedimgs );
    delete faux_flashes;
    return results;
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

    // Get the Wire Info
    TTree* fWireTree = (TTree*)fGeoFile.Get( "imagedivider/wireInfo" );
    int wireID;
    int planeID;
    float start[3];
    float end[3];
    fWireTree->SetBranchAddress( "wireID", &wireID );
    fWireTree->SetBranchAddress( "plane",  &planeID );
    fWireTree->SetBranchAddress( "start", start );
    fWireTree->SetBranchAddress( "end", end );
      
    int nentries = fWireTree->GetEntries();
    for ( int ientry=0; ientry<nentries; ientry++ ) {
      fWireTree->GetEntry(ientry);
      if ( m_WireData.find( planeID )==m_WireData.end() ) {
	// cannot find instance of wire data for plane id. make one.
	m_WireData[planeID] = larcv::pmtweights::WireData( planeID );
      }
      // start is always the one with the lowest y-value
      if ( start[1]>end[1] ) {
	// swap start end end
	for (int i=0; i<3; i++) {
	  float temp = start[i];
	  start[i] = end[i];
	  end[i] = temp;
	}
      }
      m_WireData[planeID].addWire( wireID, start, end );
    }
    
    fGeoFile.Close();

  }
  
  bool FlashMuonTaggerAlgo::findClusterEnds( const dbscan::dbscanOutput& clout, const dbscan::dbPoints& winpoints, 
					     const int clusterid, const int row_target, const int plane, 
					     const larcv::ImageMeta& meta,
					     BoundaryEndPt& endpt, larcv::Image2D& markedimg  ) {
    
    // find extremities of cluster (in time)
    int tmax = -1;
    int tmin = -1;
    int wmax = -1;
    int wmin = -1;

    for ( int ichit=0; ichit<clout.clusters.at(clusterid).size(); ichit++ ) {
      int hitidx = clout.clusters.at(clusterid).at(ichit);
      int x_ = (int)winpoints.at(hitidx).at(0)+0.1;
      int y_ = (int)winpoints.at(hitidx).at(1)+0.1;
      //markedimg.set_pixel( y_, x_, 10.0 );
      if ( tmax==-1 || y_>tmax ) { tmax = y_; wmax = x_; };
      if ( tmin==-1 || y_<tmin ) { tmin = y_; wmin = x_; };
    }
    if ( fConfig.verbosity<2 ) {
      std::cout << "end points: max (r,c)=(" << tmax << ", " << wmax << ")"
		<< " tw=(" << meta.pos_y(tmax) << "," << meta.pos_x(wmax) << ")"
		<< "  min (r,c)=(" << tmin << "," << wmin << ")" 
		<< " tw=(" << meta.pos_y(tmin) << "," << meta.pos_x(wmin) << ")"
		<< "  versus: query_row=" << row_target
		<< std::endl;
    }

    // is this a marked flash-end?
    // if extrema matches the annode flash hypothesis row. mark as interesting (Score 200)
    bool success = false;
    bool tmin_isend = false;
    bool tmax_isend = false;
    if ( abs(tmin-row_target)<=fConfig.endpoint_time_neighborhood.at(plane) ) {
      if ( fConfig.verbosity<1 ) std::cout << "MATCHED MIN END" << std::endl;
      markedimg.set_pixel( tmin, wmin, 200.0 );
      tmin_isend = true;
    }
    if ( abs(tmax-row_target)<=fConfig.endpoint_time_neighborhood.at(plane) ) { 
      if ( fConfig.verbosity<1 ) std::cout << "MATCHED MAX_END" << std::endl;
      markedimg.set_pixel( tmax, wmax, 200.0 );
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
    
    if ( success ) {
      // mark good cluster, so we don't use it again
      for ( int ichit=0; ichit<clout.clusters.at(clusterid).size(); ichit++ ) {
	int hitidx = clout.clusters.at(clusterid).at(ichit);
	int x_ = (int)winpoints.at(hitidx).at(0)+0.1;
	int y_ = (int)winpoints.at(hitidx).at(1)+0.1;
	if ( markedimg.pixel( y_, x_ )<100 )
	  markedimg.set_pixel( y_, x_, 10.0 );      
      }	      
    }
    
    return success;
  }

  float FlashMuonTaggerAlgo::calculateIntersectionTriangle( const std::vector<float>& p0, const std::vector<float>& p1, const std::vector<float>& p2 ) {

    float d01[2] = { p1[0]-p0[0], p1[1]-p0[1] }; // line segment from p0 -> p1
    float d20[2] = { p0[0]-p2[0], p0[1]-p2[1] }; // line segment from p2 -> p0
    float base = sqrt( d01[0]*d01[0] + d01[1]*d01[1] ); // length of d01, acting as our base
    float vp[2] = { d01[1], -d01[0] }; // direction perpendicular to d01 ( d01.vp = 0.0 )
    float height = 0.0;
    // project d20 onto vp to get distance from d01 line to p2
    for (int i=0; i<2; i++) height += d20[i]*vp[i]; 
    height /= base;
    
    return 0.5*std::fabs(height*base); // return area of triangle
  }

  void FlashMuonTaggerAlgo::lineSegmentIntersection2D( const std::vector< std::vector<float> >& ls1, const std::vector< std::vector<float> >& ls2, 
						       std::vector<float>& insec, int& crosses ) {
    
    //std::cout << "[lineSegmentIntersection2D] begin" << std::endl;
    //std::cout << "  testing: ls1=(" << ls1[0][0] << "," << ls1[0][1] << ") -> (" << ls1[1][0] << ","<< ls1[1][1] <<")" << std::endl;
    //std::cout << "  testing: ls2=(" << ls2[0][0] << "," << ls2[0][1] << ") -> (" << ls2[1][0] << ","<< ls2[1][1] <<")" << std::endl;
    insec.resize(2,0.0);
    float Y1 = ls1[1][1] - ls1[0][1];
    float X1 = ls1[0][0] - ls1[1][0];
    float C1 = Y1*ls1[0][0] + X1*ls1[0][1];

    float Y2 = ls2[1][1] - ls2[0][1];
    float X2 = ls2[0][0] - ls2[1][0];
    float C2 = Y2*ls2[0][0] + X2*ls2[0][1];

    float det = Y1*X2 - Y2*X1;
    if ( det==0 ) { 
      // parallel
      crosses = 0;
      //std::cout << "[lineSegmentIntersection2D] end.  parallel lines" << std::endl;
      return;
    }

    insec[0] = (X2*C1 - X1*C2)/det;
    insec[1] = (Y1*C2 - Y2*C1)/det;

    // check if interesction within line segments
    // padding needed for y-wire which is vertical
    crosses = 1;
    for (int i=0; i<2; i++) {
      if ( std::min( ls1[0][i]-0.15, ls1[1][i]-0.15 ) > insec[i] || std::max( ls1[0][i]+0.15, ls1[1][i]+0.15 )<insec[i] )
	crosses = 0;
      if ( std::min( ls2[0][i]-0.15, ls2[1][i]-0.15 ) > insec[i] || std::max( ls2[0][i]+0.15, ls2[1][i]+0.15 )<insec[i] )
	crosses = 0;
      if ( crosses==0 )
	break;
    }
    
    //std::cout << "  crosses=" << crosses << ": intersection=(" << insec[0] << "," << insec[1] << ")" << std::endl;
    //std::cout << "[lineSegmentIntersection2D] end." << std::endl;
    return;
  }

  void FlashMuonTaggerAlgo::findWireIntersections( const std::vector< std::vector<int> >& wirelists, 
						   const std::vector< std::vector<float> >& valid_range,
						   std::vector< std::vector<int> >& intersections3plane,
						   std::vector< std::vector<float> >& vertex3plane,
						   std::vector<float>& areas3plane,
						   std::vector< std::vector<int> >& intersections2plane, std::vector< std::vector<float> >& vertex2plane ) {
    // takes in lists of wire id numbers
    // returns intersection combinations, the area of the triangle they make (if 3D), 2 plane intersections
    
    // we are going to go ahead and assume 3 planes for now. we can modify for generic number of planes later (never)
    
    // was hoping that 3 plane matches the most common situation. looks like we get that if we're lucky

    // we're going to loop through all combination of wires. This can get bad fast. N^p problem...
    // expecting about O(10) endpoints per flash.  worst case, looking at 1000 combinations.
    // but we use the validity ranges to start removing intersections we don't care about

    const int nplanes = wirelists.size();
    std::set< std::vector<int> > checked2plane;
    std::set< std::vector<int> > checked3plane;

    for (int p0=0;p0<nplanes;p0++) {
      // loop through first wire plane
      for (int idx0=0; idx0<wirelists.at(p0).size(); idx0++) {
      
	// get the first wire
	int wid0 = wirelists.at(p0).at(idx0);
	const std::vector<float>& start0 = m_WireData[p0].wireStart.find(wid0)->second;
	const std::vector<float>& end0   = m_WireData[p0].wireEnd.find(wid0)->second;
	//const std::vector<float>& dir0   = m_WireData[p0].wireDir.find(wid0)->second;
	
	// go to the other planes and check the wires there
	for (int p1=p0+1; p1<nplanes; p1++) {
	  // get wire on this plane
	  for (int idx1=0; idx1<wirelists.at(p1).size(); idx1++) {
	    int wid1 = wirelists.at(p1).at(idx1);

	    std::vector< int > combo2d(3,-1);
	    combo2d[p0] = wid0;
	    combo2d[p1] = wid1;
	    if ( checked2plane.find( combo2d )==checked2plane.end() )
	      checked2plane.insert( combo2d );
	    else {
	      //std::cout << "  .. already checked: (" << combo2d[0] << "," << combo2d[1] << "," << combo2d[2] << ")" << std::endl;
	      continue;
	    }

	    const std::vector<float>& start1 = m_WireData[p1].wireStart.find(wid1)->second;
	    const std::vector<float>& end1   = m_WireData[p1].wireEnd.find(wid1)->second;
	    //const std::vector<float>& dir1   = m_WireData[p1].wireDir.find(wid1)->second;
	    
	    // change line end points from 3D to 2D (x,y,z) -> (z,y)
	    std::vector< std::vector<float> > ls0(2); // line segment 1
	    ls0[0].push_back( start0[2] );
	    ls0[0].push_back( start0[1] );
	    ls0[1].push_back( end0[2] );
	    ls0[1].push_back( end0[1] );
	    std::vector< std::vector<float> > ls1(2); // line segment 2
	    ls1[0].push_back( start1[2] );
	    ls1[0].push_back( start1[1] );
	    ls1[1].push_back( end1[2] );
	    ls1[1].push_back( end1[1] );
	    
	    // test intersection
	    std::vector<float> intersection01; 
	    int crosses01 = 0;
	    lineSegmentIntersection2D( ls0, ls1, intersection01, crosses01 );
	    if ( !crosses01 ) {
	      //std::cout << "  (" << p0 << "," << wid0 << ") and (" << p1 << "," << wid1 << ") does not cross" << std::endl;
	      continue; // move on if doesn't cross
	    }
	    bool valid = true;
	    for (int i=0; i<2; i++) {
	      if ( intersection01[i]<valid_range[i][0] || intersection01[i]>valid_range[i][1] ) {
		valid = false;
		break;
	      }
	    }
	    if ( !valid ) {
	      //std::cout << "  (" << p0 << "," << wid0 << ") and (" << p1 << "," << wid1 << ") crosses but not valid. "
	      //	<< " intersection(z,y)=(" << intersection01[0] << "," << intersection01[1] << ")"
	      //	<< std::endl;
	      continue; // not a valid intersection
	    }

	    // we got a 2plane intersection
	    std::vector<int> p2intersect(3,-1);
	    p2intersect[p0] = wid0;
	    p2intersect[p1] = wid1;
	    intersections2plane.emplace_back( p2intersect );
	    vertex2plane.push_back( intersection01 );

	    // we try for the 3plane intersection
	    int p2 = 2;
	    if ( p0==0 && p1==1 ) p2 = 2;
	    else if ( p0==0 && p1==2 ) p2 = 1;
	    else if ( p0==1 && p1==2 ) p2 = 0;
	    else
	      continue;


	    for (int idx2=0; idx2<(int)wirelists.at(p2).size(); idx2++) {
	      int wid2 = wirelists.at(p2).at(idx2);

	      std::vector< int > combo3d(3,-1);
	      combo3d[p0] = wid0;
	      combo3d[p1] = wid1;
	      combo3d[p2] = wid2;
	      if ( checked3plane.find( combo3d )==checked3plane.end() )
		checked3plane.insert( combo3d );
	      else {
		//std::cout << "    ... already checked: (" << combo3d[0] << ", " << combo3d[1] << ", " << combo3d[2] << ")" << std::endl;
		continue;
	      }

	      const std::vector<float>& start2 = m_WireData[p2].wireStart.find(wid2)->second;
	      const std::vector<float>& end2   = m_WireData[p2].wireEnd.find(wid2)->second;
	      //const std::vector<float>& dir2   = m_WireData[p2].wireDir.find(wid2)->second;

	      std::vector< std::vector<float> > ls2(2); // line segment 2
	      ls2[0].push_back( start2[2] );
	      ls2[0].push_back( start2[1] );
	      ls2[1].push_back( end2[2] );
	      ls2[1].push_back( end2[1] );

	      std::vector<float> intersection02;
	      int crosses02 = 0;
	      lineSegmentIntersection2D( ls0, ls2, intersection02, crosses02 );

	      std::vector<float> intersection12;
	      int crosses12 = 0;
	      lineSegmentIntersection2D( ls1, ls2, intersection12, crosses12 );
	      
	      if ( !crosses02 || !crosses12 )  {
		//std::cout << "  3-plane check: one combination does not cross, crosses02=" << crosses02 << " crosses12=" << crosses12 << std::endl;
		continue;
	      }

	      bool valid2 = true;
	      for (int i=0; i<2; i++) {
		if ( intersection02[i]<valid_range[i][0] || intersection02[i]>valid_range[i][1] ) {
		  //std::cout << "  3-plane check: intersection02 not valid" << std::endl;
		  valid = false;
		  break;
		}
		if ( intersection12[i]<valid_range[i][0] || intersection12[i]>valid_range[i][1] ) {
		  //std::cout << "  3-plane check: intersection12 not valid" << std::endl;
		  valid = false;
		  break;
		}
	      }
	      if ( !valid2 )
		continue; // not a valid intersection
	      
	      // got a 3 plane intersection!
	      std::vector<int> p3intersect(3,-1);
	      p3intersect[p0] = wid0;
	      p3intersect[p1] = wid1;
	      p3intersect[p2] = wid2;
	      intersections3plane.emplace_back( p3intersect );
	      // get score for 3 plane intersections
	      float area = calculateIntersectionTriangle( intersection01, intersection02, intersection12 );
	      areas3plane.push_back(area);
	      // get the 3 plane vertex
	      std::vector<float> vert3(2,0.0);
	      for (int i=0; i<2; i++) 
		vert3[i] = (intersection01[i]+intersection02[i]+intersection12[i] )/3.0;
	      vertex3plane.emplace_back( vert3 );
	    }//end of loop over p2 wires
	  }//end of loop over p1 wires
	}//end of loop over p1 planes
      }//end of loop over p0 wires
    }//end of loop over p0 planes
  }//end of findWireIntersections(...)

}

