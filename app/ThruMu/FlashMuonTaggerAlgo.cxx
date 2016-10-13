#include "FlashMuonTaggerAlgo.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

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
	if ( fSearchMode==kAnode ) {
	  tick_target = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
	}
	else if ( fSearchMode==kCathode ) {
	  tick_target = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
	  tick_target += fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick;
	}
	else if ( fSearchMode==kOutOfImage ) {
	  tick_target = opflash.Time(); // dummy opflash gives first or last tick of image
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
	  std::cout << "[Opflash search]" << std::endl;
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
    
  }//end of marked output


  bool FlashMuonTaggerAlgo::flashMatchTrackEnds( const std::vector< larlite::event_opflash* >& opflashsets, const std::vector<larcv::Image2D>& tpc_imgs,
						 std::vector< std::vector< BoundaryEndPt > >& trackendpts, std::vector< larcv::Image2D >& markedimgs ) {

    const int nplanes = tpc_imgs.size();
    
    if ( nplanes==0 )
      return false;
    
    // get a meta
    const larcv::ImageMeta& meta = tpc_imgs.at(0).meta();
    
    for ( auto& ptr_event_flash : opflashsets ) {
      for ( auto& opflash : *ptr_event_flash ) {
		
	float tick_target = 0;
	if ( fSearchMode==kAnode ) {
	  tick_target = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
	}
	else if ( fSearchMode==kCathode ) {
	  tick_target = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
	  tick_target += fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick;
	}
	else if ( fSearchMode==kOutOfImage ) {
	  tick_target = opflash.Time(); // dummy opflash gives first or last tick of image
	}
	else {
	  std::cout << "[ERROR] wrong search mode" << std::endl;
	  return false;
	}

	// check if we ned to search for this opflash
	if ( tick_target<meta.min_y() || tick_target>=meta.max_y() )
	  continue;

	// get a z-position range
	float qtot = 0;
	float z_weighted = 0.;
	for (int ipmt=0; ipmt<32; ipmt++) {
	  z_weighted += opflash.PE( ipmt )*pmtpos[ipmt][2];
	  qtot += opflash.PE( ipmt );
	}
	if ( qtot>0 ) {
	  z_weighted /= qtot;
	}
	float z_range[2] = { (float)(z_weighted-100.0), (float)(z_weighted+100.0) };

	int row_target = meta.row( tick_target );

	if ( true || fConfig.verbosity<1 ) {
	  std::cout << "============================================================================================" << std::endl;
	  std::cout << "[Opflash search]" << std::endl;
	  std::cout << "  opflash: "
		    << " tick_target=" << tick_target
		    << " row_target=" << row_target
		    << " qtot= " << qtot
		    << " z_range=[" << z_range[0] << "," << z_range[1] << "] "
		    << " drift_t=" << fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick << " ticks"
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
	
	// now that we have end points for each cluster, we look for charge deposition consistent with the flash

	// so, let's look at plane 2 first as it's dead simple to check for consistency. 
	// If there are any hits that match in z, we will make endpoint combinations with wires from the other planes that intersect it
	std::vector<int> idx_p2_flashmatched;
	for ( int e2=0; e2<endpts[2]->size(); e2++ ) {
	  float e2_w = (int)endpts[2]->at(e2).w*meta.pixel_width();
	  float e2_z = e2_w*0.3;
	  if ( z_range[0]<=e2_z && e2_z<=z_range[1] ) {
	    idx_p2_flashmatched.push_back(e2_w);
	    std::cout << "    * flash matched Y-wire: " << e2_w << " z=" << e2_z << std::endl;
	  }
	}

	// check if the U,V wires intersect with any of the y-plane matches
	std::vector<int> idx_flashmatched[2]; //< container for wires that have matched

	// we go after 3-plane intersections first. we'll do 2-plane later if we must
	std::vector< std::vector<int> > idx_3plane_intersections; // list of (u,v,y) triples that intersect
	for (int e2=0; e2<idx_p2_flashmatched.size(); e2++) {
	  // start with a y-wire and get the z-position
	  int y_wireid = idx_p2_flashmatched.at(e2);
	  const std::vector<float>& ystart = m_WireData[2].wireStart.at( y_wireid );

	  // scan for intersections in the other wire planes
	  for (int p=0; p<2; p++) {
	    const larcv::pmtweights::WireData& wiredata = m_WireData[p];
	    for (int e0=0; e0<endpts[p]->size(); e0++) {
	      float e0_w = endpts[p]->at(e0).w;
	      int wid = (int)e0_w*meta.pixel_width();
	      const std::vector<float>& start = wiredata.wireStart.find(wid)->second;
	      const std::vector<float>& end   = wiredata.wireEnd.find(wid)->second;
	      const std::vector<float>& dir0  = wiredata.wireDir.find(wid)->second;
	      
	      //const std::vector<float>& ydir   = m_WireData[2].wireDir.at( y_wireid );
	    
	      // intersection with the y-wire
	      float dz = (ystart[2]-start[2])/dir0[2];
	      float dy = dir0[1]*dz;
	      bool intersectsy = false;
	      float ypoint = start[1]+dy;
	      if ( start[1]<start[1]+dy && start[1]+dy < end[1] ) {
		idx_flashmatched[p].push_back( wid );
		std::cout << "    * flash matched plane=" << p << ": intersection (z,y)=(" << ystart[2] << ", " << start[1]+dy << ")" << std::endl;
		intersectsy = true;
	      }

	      // move 
	      if ( !intersectsy )
		continue;
	    }
	  }
	}//end of loop over all planes by the Y-plane
	
      }//end of opflashes loop
    }//end of opflashsets    
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
      success =true;
      markedimg.set_pixel( tmax, wmax, 200.0 );
      tmin_isend = true;
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

}

