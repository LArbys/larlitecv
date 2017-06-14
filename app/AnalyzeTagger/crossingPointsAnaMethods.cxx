#include "crossingPointsAnaMethods.h"

#include "TTree.h"
#include "TLorentzVector.h"

// larcv
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/ImageMeta.h"
#include "SCE/SpaceChargeMicroBooNE.h"
#include "UBWireTool/UBWireTool.h"

// larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/trigger.h"

// larlitecv
#include "TaggerTypes/dwall.h"

namespace larlitecv {

  void CrossingPointAnaData_t::bindToTree( TTree* tree ) {
    tree->Branch( "tot_proposed_crossingpoints", &tot_proposed_crossingpoints, "tot_proposed_crossingpoints/I" );
    tree->Branch( "tot_true_crossingpoints", &tot_true_crossingpoints, "tot_true_crossingpoints/I" );
    tree->Branch( "tot_flashmatched_true_crossingpoints", &tot_flashmatched_true_crossingpoints, "tot_flashmatched_true_crossingpoints/I" );
    tree->Branch( "tot_matched_crossingpoints", &tot_matched_crossingpoints, "tot_matched_crossingpoints/I" );    
    tree->Branch( "proposed_crossingpoints", proposed_crossingpoints, "proposed_crossingpoints[7]/I" );
    tree->Branch( "true_crossingpoints",     true_crossingpoints, "true_crossingpoints[6]/I" );
    tree->Branch( "flashmatched_true_crossingpoints", flashmatched_true_crossingpoints, "flashmatched_true_crossingpoints[6]/I" );
    tree->Branch( "matched_crossingpoints",  matched_crossingpoints, "matched_crossingpoints[6]/I" );
    tree->Branch( "true_intime_thrumu",      &true_intime_thrumu,     "true_intime_thrumu/I");
    tree->Branch( "true_intime_stopmu",      &true_intime_stopmu,     "true_intime_stopmu/I");
  }

  void analyzeCrossingPoints( CrossingPointAnaData_t& data, const larcv::ImageMeta& meta,
			      const larcv::EventPixel2D* ev_spacepoints[], const std::vector<larlite::event_opflash*>& opflash_v,
			      const larlite::trigger* ev_trigger, const larlite::event_mctrack* ev_mctrack ) {
			      
    data.clear();

    if ( ev_mctrack && ev_trigger ) {
      analyzeCrossingMCTracks( data, meta, ev_trigger, ev_mctrack, opflash_v, true );
    }
  }
  
  void analyzeCrossingMCTracks( CrossingPointAnaData_t& data, const larcv::ImageMeta& meta,
				const larlite::trigger* ev_trigger, const larlite::event_mctrack* ev_mctrack, const std::vector<larlite::event_opflash*>& opflash_v,
				bool printFlashEnds ) {
    // loop over MC tracks, get end points of muons
    //int intime_cosmics = 0;

    larlitecv::SpaceChargeMicroBooNE sce;
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;    
    
    for ( auto const& track : *ev_mctrack ) {
      if ( std::abs(track.PdgCode())!=13  ) continue; // skip muons for now
      if ( track.size()==0 ) continue;
      const TLorentzVector& track_start = track.front().Position();
      std::vector<float> fstart(3,0);
      fstart[0] = track_start.X();
      fstart[1] = track_start.Y();
      fstart[2] = track_start.Z();
      std::vector<float> fend(3,0);
      fend[0]   = track.back().Position().X();
      fend[1]   = track.back().Position().Y();
      fend[2]   = track.back().Position().Z();
      
      int track_start_boundary = 0;
      float track_start_dwall = larlitecv::dwall( fstart, track_start_boundary );
      int track_end_boundary = 0;
      float track_end_dwall = larlitecv::dwall( fend, track_end_boundary );

      const larlite::mcstep& first_step = track.front();
      const larlite::mcstep& last_step  = track.back();

      // space charge corrections
      std::vector<double> start_offset = sce.GetPosOffsets( first_step.X(), first_step.Y(), first_step.Z() );
      std::vector<double> end_offset   = sce.GetPosOffsets( last_step.X(),  last_step.Y(),  last_step.Z() );

      std::vector<float> sce_start(3);
      sce_start[0] = first_step.X()-start_offset[0]+0.7;
      sce_start[1] = first_step.Y()+start_offset[1];
      sce_start[2] = first_step.Z()+start_offset[2];
      Double_t sce_start_xyz[3] = { sce_start[0], sce_start[1], sce_start[2] };

      std::vector<float> sce_end(3);
      sce_end[0] = last_step.X()-end_offset[0]+0.7;
      sce_end[1] = last_step.Y()+end_offset[1];
      sce_end[2] = last_step.Z()+end_offset[2];
      Double_t sce_end_xyz[3] = { sce_end[0], sce_end[1], sce_end[2] };

      std::vector< int > start_pix(4); // (row, u-plane, v-plane, y-plane)
      std::vector< int > end_pix(4); // (row, u-plane, v-plane, y-plane)
      start_pix[0] = (first_step.T()*1.0e-3 - (ev_trigger->TriggerTime()-4050.0))/0.5 + sce_start[0]/cm_per_tick + 3200;
      end_pix[0]   = (last_step.T()*1.0e-3  - (ev_trigger->TriggerTime()-4050.0))/0.5 + sce_end[0]/cm_per_tick   + 3200;
      for ( size_t p=0; p<3; p++ ) {
	start_pix[p+1] = larutil::Geometry::GetME()->WireCoordinate( sce_start_xyz, p );
	end_pix[p+1]   = larutil::Geometry::GetME()->WireCoordinate( sce_end_xyz,   p );
      }

      float orig_start_tick = (first_step.T()*1.0e-3 - (ev_trigger->TriggerTime()-4050.0))/0.5 + 3200.0;
      float orig_end_tick   = (last_step.T()*1.0e-3  - (ev_trigger->TriggerTime()-4050.0))/0.5 + 3200.0;
      float drift_ticks      = 258.0/cm_per_tick;
      bool start_intime = false;
      bool start_crosses = false;
      if ( start_pix[0]>meta.min_y() && start_pix[0]<meta.max_y() ) {
	start_intime = true;
	start_pix[0] = meta.row( start_pix[0] );
	  
	if ( track_start_dwall < 10.0 ) {

	  if ( printFlashEnds && track_start_boundary==4 )
	    std::cout << "Anode Boundary Crossing:   row=" << start_pix[0] << " tick=" << meta.pos_y(start_pix[0]) << " (orig=" << orig_start_tick << ")"
		      << " pos=(" << first_step.X() << "," << first_step.Y() << "," << first_step.Z() << ")" << std::endl;
	  else if ( printFlashEnds && track_start_boundary==5 )
	    std::cout << "Cathode Boundary Crossing: row=" << start_pix[0] << " tick=" << meta.pos_y(start_pix[0]) << " (flash=" << orig_start_tick << " orig=" << orig_start_tick+drift_ticks << ")"
		      << " pos=(" << first_step.X() << "," << first_step.Y() << "," << first_step.Z() << ")" << std::endl;

	  // flash match this crossing point
	  bool found_flash = false;	  
	  if ( track_start_boundary==4 || track_start_boundary==5 ) {
	    for ( auto const& ev_opflash : opflash_v ) {
	      for ( auto const& opflash : *ev_opflash ) {
		if ( found_flash )
		  break;		
		if ( fabs( orig_start_tick - (3200.0+opflash.Time()/0.5) )<5.0 ) {
		  found_flash = true;
		}
	      }
	      if ( found_flash )
		break;
	    }//end of ev_opflash loop
	  }
	  else
	    found_flash = true; // no flash to match for others, but want to count them in array for convenience

	  if ( found_flash ) {
	    data.tot_flashmatched_true_crossingpoints++;
	    data.flashmatched_true_crossingpoints[ track_start_boundary ]++;
	  }

	  data.start_pixels.emplace_back( std::move(start_pix) );
	  data.start_crossingpts.emplace_back( std::move(sce_start) );
	  data.start_type.push_back( track_start_boundary );
	  start_crosses = true;
	  data.tot_true_crossingpoints++;
	  data.true_crossingpoints[track_start_boundary]++;
	}
      }
      
      bool end_intime = false;
      bool end_crosses = false;
      if ( end_pix[0]>meta.min_y() && end_pix[0]<meta.max_y() ) {
	end_intime = true;
	end_pix[0]   = meta.row( end_pix[0] );
	if ( track_end_dwall < 10.0 ) {

	  if ( printFlashEnds && track_end_boundary==4 )
	    std::cout << "Anode Boundary Crossing:   row=" << end_pix[0] << " tick=" << meta.pos_y(end_pix[0]) << " (orig=" << orig_end_tick << ")"
		      << " pos=(" << last_step.X() << "," << last_step.Y() << "," << last_step.Z() << ")" << std::endl;
	  else if ( printFlashEnds && track_end_boundary==5 )
	    std::cout << "Cathode Boundary Crossing: row=" << end_pix[0] << " tick=" << meta.pos_y(end_pix[0])  << " (flash=" << orig_end_tick << " orig=" << orig_end_tick+drift_ticks << ")"
		      << " pos=(" << last_step.X() << "," << last_step.Y() << "," << last_step.Z() << ")" << std::endl;

	  // flash match this crossing point
	  bool found_flash = false;	  
	  if ( track_end_boundary==4 || track_end_boundary==5 ) {
	    for ( auto const& ev_opflash : opflash_v ) {
	      for ( auto const& opflash : *ev_opflash ) {
		if ( found_flash )
		  break;		
		if ( fabs( orig_end_tick - (3200.0+opflash.Time()/0.5) )<5.0 ) {
		  found_flash = true;
		}
	      }
	      if ( found_flash )
		break;
	    }//end of ev_opflash loop
	  }
	  else
	    found_flash = true; // no flash to match for others, but want to count them in array for convenience
	  
	  if ( found_flash ) {
	    data.tot_flashmatched_true_crossingpoints++;
	    data.flashmatched_true_crossingpoints[ track_end_boundary ]++;
	  }
	  
	  data.end_pixels.emplace_back( std::move(end_pix) );
	  data.end_crossingpts.emplace_back( std::move(sce_end) );
	  data.end_type.push_back( track_end_boundary );
	  end_crosses = true;
	  data.tot_true_crossingpoints++;
	  data.true_crossingpoints[track_end_boundary]++;
	}
      }
      
      if ( start_intime || end_intime ) {
	//intime_cosmics++;
	if ( start_intime && !start_crosses ) {
	  std::cout << "start point does not cross boundary: (" << fstart[0] << "," << fstart[1] << "," << fstart[2] << ")"
		    << " dwall=" << track_start_dwall << " intime=" << start_intime
		    << " tick=" << meta.pos_y( start_pix[0] )
		    << " row=" << start_pix[0]
		    << std::endl;
	  //throw std::runtime_error("start point does not cross boundary?");
	}
	else if ( start_crosses && end_crosses )
	  data.true_intime_thrumu++;
	else if ( start_crosses && !end_crosses )
	  data.true_intime_stopmu++;
      }	
      
    }//end of loop over mctracks
    std::cout << "number of intime thrumu: "        << data.true_intime_thrumu << std::endl;
    std::cout << "number of intime stopmu: "        << data.true_intime_stopmu << std::endl;
    std::cout << "number of true crossing points: " << data.tot_true_crossingpoints << std::endl;
  }

  void analyzeCrossingDataOnly( CrossingPointAnaData_t& data, std::vector<larcv::EventPixel2D*>& ev_spacepoints ) {
    //analyze proposed boundary points
    //std::cout << "Analyze Boundary Points" << std::endl;
    data.tot_proposed_crossingpoints = 0;
    for (int i=0; i<7; i++) {
      if ( ev_spacepoints[i]==NULL)
	throw std::runtime_error("wtf");
      //std::cout << " endtype " << spacepoint_producers[i] << ": ";
      std::cout << ev_spacepoints[i]->Pixel2DArray(0).size() << std::endl;
      data.proposed_crossingpoints[i]  += ev_spacepoints[i]->Pixel2DArray(0).size();
      data.tot_proposed_crossingpoints += ev_spacepoints[i]->Pixel2DArray(0).size();
    }
  }

  void analyzeCrossingMatches( CrossingPointAnaData_t& data, std::vector<larcv::EventPixel2D*> ev_spacepoints, const larcv::ImageMeta& meta ) {
    // MC info already stored in data

    std::vector<bool> matched_startpoint( data.start_pixels.size(), false );
    std::vector<bool> matched_endpoint(   data.end_pixels.size(),   false );

    std::vector< std::vector<int> >* p_crossing_pixel_v[2] = { &data.start_pixels, &data.end_pixels }; // truth pixels for crossing points
    std::vector< std::vector<float> >* p_crossingpts_v[2]  = { &data.start_crossingpts, &data.end_crossingpts }; // truth 3D positions
    std::vector<bool>* p_matched_v[2] = { &matched_startpoint, &matched_endpoint };
    
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;

    for ( int v=0; v<2; v++ ) {
      for ( int ipix=0; ipix<(int)p_crossing_pixel_v[v]->size(); ipix++ ) {

	// we need to get the 3D position to compare against
	const std::vector<int>&  pixinfo     = p_crossing_pixel_v[v]->at(ipix); // position in image
	const std::vector<float>& crossingpt = p_crossingpts_v[v]->at(ipix);    // position in 3D

	// use TPC position to get X
	std::vector<float> crossingpt_tpcx(3,0);
	crossingpt_tpcx[0] = (meta.pos_y( pixinfo[0] )-3200.0)*cm_per_tick;
	crossingpt_tpcx[1] = crossingpt[1];
	crossingpt_tpcx[2] = crossingpt[2];

	// scan for pixel, loop over types and pts
	bool matched = false;
	for (int i=0; i<6; i++) {
	  if ( matched )
	    break;

	  for ( int j=0; j<(int)ev_spacepoints[i]->Pixel2DArray(0).size(); j++ ) {

	    std::vector<float> intersect(2,0.0);
	    std::vector<int> wids(3,0);
	    int crossing = 0;
	    double triangle_area = 0.0;
	    for (int p=0; p<3; p++) {
	      wids[p] = ev_spacepoints[i]->Pixel2DArray(p).at(j).X();
	    }
	    

	    larcv::UBWireTool::wireIntersection( wids, intersect, triangle_area, crossing );

	    if ( crossing==0 )
	      continue;
	    
	    float x = ( meta.pos_y( ev_spacepoints[i]->Pixel2DArray(0).at(j).Y() ) - 3200.0 )*cm_per_tick;

	    std::vector<float> spacepoints(3);
	    spacepoints[0] = x;
	    spacepoints[1] = intersect[1];
	    spacepoints[2] = intersect[0];

	    float dist = 0;
	    for (int d=0; d<3; d++) {
	      dist += (spacepoints[d]-crossingpt_tpcx[d])*(spacepoints[d]-crossingpt_tpcx[d]);
	    }
	    dist = sqrt(dist);
	    //std::cout << "true[" << v << "," << ipix << "] vs. proposed[" << i << "," << j << "] dist=" << dist << std::endl;
	    if (dist<20.0) {
	      matched = true;
	    }

	    if ( matched )
	      break;
	  }// end of loop over tagged points of type i
	}//end of boundary point types

	p_matched_v[v]->at(ipix) = matched;
	
      }
    }//end of v loop

    for (size_t i=0; i<data.start_type.size(); i++) {
      if ( matched_startpoint[i] ) {
	data.matched_crossingpoints[ data.start_type[i] ]++;
	data.tot_matched_crossingpoints++;
      }
    }

    for (size_t i=0; i<data.end_type.size(); i++) {
      if ( matched_endpoint[i] ) {
	data.matched_crossingpoints[ data.end_type[i] ]++;
	data.tot_matched_crossingpoints++;
      }
    }
    std::cout << "Proposed Crossing Points: " << data.tot_proposed_crossingpoints << std::endl;
    std::cout << "True Crossing Points: "     << data.tot_true_crossingpoints << std::endl;    
    std::cout << "Matched Crossing Points: "  << data.tot_matched_crossingpoints << std::endl;
    
  }
  
}
