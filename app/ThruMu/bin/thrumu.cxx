#include <iostream>
#include <cmath>

// config/storage
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite data
#include "Base/DataFormatConstants.h"
#include "DataFormat/opflash.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "ANN/ANNAlgo.h"
#include "dbscan/DBSCANAlgo.h"

// larelite
#include "ThruMu/BoundaryMuonTaggerAlgo.h"

// algos



int main( int nargs, char** argv ) {
  
  std::cout << "[BOUNDARY MUON TAGGER]" << std::endl;
  
  // configuration
  larcv::PSet cfg = larcv::CreatePSetFromFile( "bmt.cfg" );
  larcv::PSet bmt = cfg.get<larcv::PSet>("BoundaryMuonTagger");
  larcv::PSet sidetagger_pset  = bmt.get<larcv::PSet>("BMTSideTagger");
  larcv::PSet flashtagger_pset = bmt.get<larcv::PSet>("BMTFlashTagger");

  std::string larcv_image_producer = bmt.get<std::string>("InputLArCVImages");
  
  // Configure Data coordinator
  larlitecv::DataCoordinator dataco;

  // larlite
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_mcinfo_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_wire_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_opdigit_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_opreco_0000.root", "larlite" );

  // larcv
  dataco.add_inputfile( "data/data_samples/v05/spoon/larcv/spoon_larcv_out_0000.root", "larcv" );

  // configure
  dataco.configure( "bmt.cfg", "StorageManager", "IOManager", "BoundaryMuonTagger" );
  
  // initialize
  dataco.initialize();


  // Configure Algorithms
  // side-tagger
  larlitecv::BoundaryMuonTaggerAlgo sidetagger;
  larlitecv::ConfigBoundaryMuonTaggerAlgo sidetagger_cfg;
  sidetagger_cfg.neighborhoods = sidetagger_pset.get< std::vector<int> >("Neighborhoods");
  sidetagger_cfg.thresholds    = sidetagger_pset.get< std::vector<float> >( "Thresholds" );
  sidetagger.configure(sidetagger_cfg);

  // flash-tagger
  std::vector<std::string> opflash_producers = flashtagger_pset.get< std::vector<std::string> >( "OpflashProducers" );
  std::vector< int > time_neighborhood       = flashtagger_pset.get< std::vector<int> >( "TimeNeighboorHood" ); // one for each plane
  std::vector< int > wire_neighborhood       = flashtagger_pset.get< std::vector<int> >( "WireNeighboorHood" ); // one for each plane
  std::vector< float > charge_threshold      = flashtagger_pset.get< std::vector<float> >( "ChargeThreshold" ); // one for each plane


  // Start Event Loop
  //int nentries = dataco.get_nentries("larcv");
  //int nentries = 5;
  int nentries = 1;
  for (int ientry=0; ientry<nentries; ientry++) {
    
    dataco.goto_entry(ientry,"larcv");

    // get images (from larcv)
    std::cout << "get data" << std::endl;
    larcv::EventImage2D* event_imgs    = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, larcv_image_producer );

    // ------------------------------------------------------------------------------------------//
    // SIDE TAGGER //
    // first we run the boundary muon tagger algo that matches top/side/upstream/downstream
    std::cout << "search for boundary pixels" << std::endl;
    std::vector< larcv::Image2D > outhits;
    sidetagger.searchforboundarypixels( event_imgs->Image2DArray(), outhits );

    // ------------------------------------------------------------------------------------------//
    // FLASH TAGGER //

    // create storage for new images
    std::vector< larcv::Image2D > flashtagger_hits;

    // new image for flash hits
    std::vector<larcv::Image2D> stage1_annode_hits;  // all in-time hits (non-clustered, non-edged)
    std::vector<larcv::Image2D> stage1_cathode_hits; // all in-time hits (non-clustered, non-edged)
    for ( auto &tpc_img : event_imgs->Image2DArray() ) {
      larcv::Image2D annode_img( tpc_img.meta() );
      larcv::Image2D cathode_img( tpc_img.meta() );
      annode_img.paint(0.0);
      cathode_img.paint(0.0);
      stage1_annode_hits.emplace_back( annode_img );
      stage1_cathode_hits.emplace_back( cathode_img );
    }

    std::vector<larcv::Image2D> annode_hits;
    std::vector<larcv::Image2D> cathode_hits;
    
    // loop through flash producers
    for ( auto &flashproducer : flashtagger_pset.get< std::vector<std::string> >( "OpflashProducers" ) ) {
      larlite::event_opflash* opdata = (larlite::event_opflash*)dataco.get_larlite_data(larlite::data::kOpFlash, flashproducer );
      std::cout << "search for flash hits from " << flashproducer << ": " << opdata->size() << " flashes" << std::endl;

      // loop through flashes
      for ( auto &opflash : *opdata ) {
	//bool hitfound[3] = {false,false,false};

	std::cout << "[opflash]" << std::endl;

	// go through each plane and tag pixels that in time with the flash (and the drift)
	for ( auto &tpc_img : event_imgs->Image2DArray() ) {
	  
	  const larcv::ImageMeta& meta = tpc_img.meta();
	  int plane = meta.plane();
      
	  float tick_annode = 3200.0 + opflash.Time()/0.5;  // (trigger tick) + (flash time from trigger)/(us per tick)
	  int   row_annode;
	  if ( tick_annode<2400 || tick_annode>=meta.max_y() )
	    row_annode = -1;
	  else
	    row_annode = meta.row( tick_annode );

	  float dtick_drift = 250.0/0.111/0.5; // cm / (cm/usec) / (usec/tick)
	  float tick_cathode = 3200.0 + opflash.Time()/0.5 + dtick_drift;  // (trigger tick) + (flash time from trigger)/(us per tick)
	  int   row_cathode;
	  if ( tick_cathode<2400 || tick_cathode>=meta.max_y() )
	    row_cathode = -1;
	  else
	    row_cathode = meta.row( tick_cathode );

	  if ( plane==0 )
	    std::cout << "  opflash: p= " << plane 
		      << " tick_annode=" << tick_annode 
		      << " row_annode=" << row_annode
		      << " tick_cathod=" << tick_cathode
		      << " row_cathode=" << row_cathode << std::endl;
	  
	  // search strategy:
	  // (1) we get a row to search, 
	  // (2) we scan across the wires until we see a pixel with charge above threshold
	  // (3) we cut out a window in time and charge
	  // (4) we scan across it, gathering above threshold hits into a point list.
	  // (5) we cluster in the point list using dbscan
	  // (6) the most extreme pixels in time are marked as endpoints
	  // (7) is the endpoint close to the flash? if so mark it.
	  // (8) in stage 2 image, we mark (10) if explored, 200 if end-tag
	  
	  // search in time and wire neighborhoods
	  // ANNODE
	  if ( row_annode!=-1 ) {
	    for ( int dt=-time_neighborhood.at(plane); dt<=time_neighborhood.at(plane); dt++ ) {
	      int r_annode  = row_annode+dt;
	      
	      // scan wires
	      for (int iwire=0; iwire<meta.cols(); iwire++) {

		// is it a good pixel to check around?
		if ( r_annode>=0 && r_annode<meta.rows()  // within range
		     && stage1_annode_hits.at(plane).pixel( r_annode, iwire) <= 1.0 // not visited
		     && tpc_img.pixel( r_annode, iwire )>charge_threshold[plane] ) {  // above threshold
		  
		  // new, unexplored region!
		  // we define a window around this point: (r_annode, c)
		  float t1 = meta.pos_y(r_annode);
		  float w  = meta.pos_x( iwire );
		  dbscan::dbPoints winpoints;
		  int centeridx = -1;
		  for (int dwire=-10; dwire<=10; dwire++) {
		    for (int dtwin=-20; dtwin<=20; dtwin++) {
		      int r = r_annode+dtwin;
		      int c = iwire+dwire;
		      // check validity
		      if ( r>=0 && r<meta.rows() && c>=0 && c<meta.cols() 
			   && tpc_img.pixel(r,c)>charge_threshold[plane] ) {
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
		  std::cout << "exploring around: (c,r)=" << iwire << ", " << r_annode << " (w,t)=(" << w << "," << t1 << ")  npoints=" << winpoints.size() << std::endl;

		  // cluster
		  dbscan::DBSCANAlgo dbalgo;
		  dbscan::dbscanOutput clout = dbalgo.scan( winpoints, 3, 3.0, false, 0.0 );

		  // which cluster is the center point in?
		  int connected_cluster = clout.clusterid.at(centeridx);
		  std::cout << "  connected clusterid=" << connected_cluster  << std::endl;
		  for (int ic=0; ic<clout.clusters.size(); ic++) {
		    std::cout << "    clusterid=" << ic << ", size=" << clout.clusters.at(ic).size() << std::endl;
		  }
		  
		  if ( clout.clusters.at(connected_cluster).size()<5 )
		    continue;

		  std::cout << "valid connected cluster: " << connected_cluster << " cluster size=" << clout.clusters.at(connected_cluster).size() << std::endl;

		  // find extremities of cluster
		  int tmax = -1;
		  int tmin = -1;
		  int wmax = -1;
		  int wmin = -1;
		  for ( int ichit=0; ichit<clout.clusters.at(connected_cluster).size(); ichit++ ) {
		    int hitidx = clout.clusters.at(connected_cluster).at(ichit);
		    int x_ = (int)winpoints.at(hitidx).at(0)+0.1;
		    int y_ = (int)winpoints.at(hitidx).at(1)+0.1;
		    stage1_annode_hits.at(plane).set_pixel( y_, x_, 10.0 );
		    if ( tmax==-1 || y_>tmax ) { tmax = y_; wmax = x_; };
		    if ( tmin==-1 || y_<tmin ) { tmin = y_; wmin = x_; };
		  }
		  std::cout << "end points: max (r,c)=(" << tmax << ", " << wmax << ")"
			    << " tw=(" << meta.pos_y(tmax) << "," << meta.pos_x(wmax) << ")"
			    << "  min (r,c)=(" << tmin << "," << wmin << ")" 
			    << " tw=(" << meta.pos_y(tmin) << "," << meta.pos_x(wmin) << ")"
			    << "  versus: r_annode=" << r_annode
			    << std::endl;
		  
		  // extrema matches the annode flash hypothesis row. mark as interesting (Score 200)
		  if ( abs(tmin-r_annode)<=time_neighborhood.at(plane) ) {
		    std::cout << "MATCHED MIN END" << std::endl;
		    stage1_annode_hits.at(plane).set_pixel( tmin, wmin, 200.0 );
		  }
		  else if ( abs(tmax-r_annode)<=time_neighborhood.at(plane) ) { 
		    std::cout << "MATCHED MAX_END" << std::endl;
		    stage1_annode_hits.at(plane).set_pixel( tmax, wmax, 200.0 );
		  }
		  
		}//if pixel valid, unexplored and above threshold
	      }//end of wire scan
	    } //loop over time-neighborhoods
	  }// if annode time tick within image
	  
	}//end of loop over images
	
      }// loop over flashes
      
    }// loop over producer
    
    // ------------------------------------------------------------------------------------------//
    // SAVE IMAGES //
    
    // from side tagger
    larcv::EventImage2D* boundary_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, sidetagger_pset.get<std::string>("OutputMatchedPixelImage") );
    boundary_imgs->Emplace( std::move(outhits) );    
    
    // flash tagger
    larcv::EventImage2D* stage1_annode_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stage1_annode" );
    larcv::EventImage2D* stage1_cathode_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stage1_cathode" );
    stage1_annode_imgs->Emplace( std::move(stage1_annode_hits) );
    stage1_cathode_imgs->Emplace( std::move(stage1_cathode_hits) );
    
    // go to tree
    dataco.save_entry();
  }
  
  dataco.finalize();

  return 0;
}
