#include "crossingPointsAnaMethods.h"

#include <cstring>
#include <stdexcept>

#include "TTree.h"
#include "TLorentzVector.h"

// larcv
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/Image2D.h"
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
#include "TaggerTypes/BoundaryMuonTaggerTypes.h"

namespace larlitecv {

  void CrossingPointAnaData_t::bindToTree( TTree* tree ) {
    tree->Branch( "tot_proposed_crossingpoints", &tot_proposed_crossingpoints, "tot_proposed_crossingpoints/I" );
    tree->Branch( "tot_true_crossingpoints", &tot_true_crossingpoints, "tot_true_crossingpoints/I" );
    tree->Branch( "tot_flashmatched_true_crossingpoints", &tot_flashmatched_true_crossingpoints, "tot_flashmatched_true_crossingpoints/I" );
    tree->Branch( "tot_matched_crossingpoints", &tot_matched_crossingpoints, "tot_matched_crossingpoints/I" );    
    tree->Branch( "proposed_crossingpoints", proposed_crossingpoints, "proposed_crossingpoints[7]/I" );
    tree->Branch( "true_crossingpoints",     true_crossingpoints, "true_crossingpoints[7]/I" );
    tree->Branch( "flashmatched_true_crossingpoints", flashmatched_true_crossingpoints, "flashmatched_true_crossingpoints[7]/I" );
    tree->Branch( "matched_crossingpoints",  matched_crossingpoints, "matched_crossingpoints[7]/I" );
    tree->Branch( "true_intime_thrumu",      &true_intime_thrumu,     "true_intime_thrumu/I");
    tree->Branch( "true_intime_stopmu",      &true_intime_stopmu,     "true_intime_stopmu/I");
  }

  void analyzeCrossingPoints( CrossingPointAnaData_t& data, const larcv::ImageMeta& meta, const std::vector<larcv::Image2D>& img_v,
			      const larcv::EventPixel2D* ev_spacepoints[], const std::vector<larlite::event_opflash*>& opflash_v,
			      const larlite::trigger* ev_trigger, const larlite::event_mctrack* ev_mctrack ) {
			      
    data.clear();

    if ( ev_mctrack && ev_trigger ) {
      analyzeCrossingMCTracks( data, meta, img_v, ev_trigger, ev_mctrack, opflash_v, true );
    }
  }
  
  void analyzeCrossingMCTracks( CrossingPointAnaData_t& data, const larcv::ImageMeta& meta, const std::vector<larcv::Image2D>& img_v,
				const larlite::trigger* ev_trigger, const larlite::event_mctrack* ev_mctrack, const std::vector<larlite::event_opflash*>& opflash_v,
				bool printFlashEnds ) {
    // loop over MC tracks, get end points of muons
    //int intime_cosmics = 0;

    larlitecv::SpaceChargeMicroBooNE sce;
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    const int endpt_radius = 10;

    data.mctrack_imgendpoint_indices.resize( ev_mctrack->size() );

    int trackindex = -1;
    for ( auto const& track : *ev_mctrack ) {
      trackindex++;
      
      if ( std::abs(track.PdgCode())!=13  ) continue; // skip muons for now
      if ( track.size()==0 ) continue; // skip garbage

      if ( printFlashEnds )
	std::cout << "[TRACK #" << trackindex << "] " << std::endl;
      

      // get the first and last steps
      const larlite::mcstep& first_step = track.front();
      const larlite::mcstep& last_step  = track.back();

      // get the position of the track ends in terms of the image tick axis
      float tick_start = getTick( first_step, ev_trigger->TriggerTime(), &sce );
      float tick_end   = getTick( last_step,  ev_trigger->TriggerTime(), &sce );

      // if the track is entirely outside the image, skip it
      if ( tick_start<meta.min_y() && tick_end<meta.min_y() ) {
	if ( printFlashEnds )
	  std::cout << "  Entire track before the image time region." << std::endl;
	continue;
      }
      if ( tick_start>meta.max_y() && tick_end>meta.max_y() ) {
	if ( printFlashEnds ) 
	  std::cout << "  Entire track after the image time region." << std::endl;
	continue;
      }
      
      // follow the track and get the first positions (detector coordinates) where the track shows up in the image
      std::vector<float> fstart = getFirstStepPosInsideImage( track, meta, ev_trigger->TriggerTime(), true, 0.15, 3.0, &sce );
      std::vector<float> fend   = getFirstStepPosInsideImage( track, meta, ev_trigger->TriggerTime(), false, 0.15, 3.0, &sce );
      
      if ( fstart.size()==0 || fend.size()==0 ) {
	if ( printFlashEnds )
	  std::cout << "  No part of track inside image." << std::endl;
	continue;
      }

      // get the boundary the track ends are close to
      int track_start_boundary = 0;
      float track_start_dwall = larlitecv::dwall( fstart, track_start_boundary );
      std::string start_crossingname = larlitecv::BoundaryEndNames( (larlitecv::BoundaryEnd_t)track_start_boundary );
	  
      int track_end_boundary = 0;
      float track_end_dwall = larlitecv::dwall( fend, track_end_boundary );
      std::string end_crossingname = larlitecv::BoundaryEndNames( (BoundaryEnd_t)track_end_boundary );      

      // space charge corrections to for fstart and fend (needed to project into the image properly)
      std::vector<double> start_offset = sce.GetPosOffsets( fstart[0], fstart[1], fstart[2] );
      std::vector<double> end_offset   = sce.GetPosOffsets( fend[0], fend[1], fend[2] );

      // get the space-charge corrected start and end points
      std::vector<float> sce_start(3);
      sce_start[0] = fstart[0]-start_offset[0]+0.7;
      sce_start[1] = fstart[1]+start_offset[1];
      sce_start[2] = fstart[2]+start_offset[2];
      Double_t sce_start_xyz[3] = { sce_start[0], sce_start[1], sce_start[2] };

      std::vector<float> sce_end(3);
      sce_end[0] = fend[0]-end_offset[0]+0.7;
      sce_end[1] = fend[1]+end_offset[1];
      sce_end[2] = fend[2]+end_offset[2];
      Double_t sce_end_xyz[3] = { sce_end[0], sce_end[1], sce_end[2] };

      // get the projected pixel location in the image
      std::vector< int > start_pix(4,0); // (row, u-plane, v-plane, y-plane)
      std::vector< int > end_pix(4,0); // (row, u-plane, v-plane, y-plane)
      start_pix[0] = (first_step.T()*1.0e-3 - (ev_trigger->TriggerTime()-4050.0))/0.5 + sce_start[0]/cm_per_tick + 3200; // tick
      end_pix[0]   = (last_step.T()*1.0e-3  - (ev_trigger->TriggerTime()-4050.0))/0.5 + sce_end[0]/cm_per_tick   + 3200; // tick
      for ( size_t p=0; p<3; p++ ) {
	start_pix[p+1] = larutil::Geometry::GetME()->WireCoordinate( sce_start_xyz, p );
	end_pix[p+1]   = larutil::Geometry::GetME()->WireCoordinate( sce_end_xyz,   p );
      }

      // does track cross the image boundary?
      // -1 does not cross
      // 0 start out -> end in
      // 1 start in -> end out
      int crosses_image_boundary = doesTrackCrossImageBoundary( track, meta, ev_trigger->TriggerTime(), &sce );      
      
      // get the time of the flash in the image clock (i.e. ticks). do not correct for time to drift to wire planes
      float orig_start_tick = (first_step.T()*1.0e-3 - (ev_trigger->TriggerTime()-4050.0))/0.5 + 3200.0;
      float orig_end_tick   = (last_step.T()*1.0e-3  - (ev_trigger->TriggerTime()-4050.0))/0.5 + 3200.0;
      float drift_ticks      = 258.0/cm_per_tick;

      // ----------------------------------------------------------------------
      // TRACK START 
      bool start_intime = false;  // start is in the image (true by construction)
      bool start_crosses = false; // crosses detector boundary
      bool start_found_flash = false; // does a flash pair with the time of the event?
      int flash_index = -1; // if a flash is paired, its index in the vector(s) are stored here
      int start_nplanes_above_thresh = 0; // number of planes where region of projected image pixel has charge above threshold

      // flash match start crossing point
      int iflashindex = -1;
      for ( auto const& ev_opflash : opflash_v ) {
	for ( auto const& opflash : *ev_opflash ) {
	  iflashindex++;
	  if ( start_found_flash )
	    break;		
	  if ( fabs( orig_start_tick - (3200.0+opflash.Time()/0.5) )<5.0 ) {
	    start_found_flash = true;
	    flash_index = iflashindex;
	  }
	}
      }

      // change from tick to row coordinate
      start_pix[0] = meta.row( start_pix[0] );
	  
      // check neighborhood around pixel (has charge around?)
      start_nplanes_above_thresh = 0;
      for (int p=0; p<3; p++) {
	int nabove_thresh = 0;
	for (int dr=-endpt_radius; dr<=endpt_radius; dr++) {
	  int srow = start_pix[0]+dr;
	  if ( srow<0 || srow>=(int)meta.rows() )
	    continue;
	  for (int dc=-endpt_radius; dc<=endpt_radius; dc++) {
	    int scol = start_pix[p+1]+dc;
	    if ( scol<0 || scol>=(int)meta.cols() )
	      continue;
	    if ( img_v[p].pixel( srow, scol )>10.0 ) {
	      nabove_thresh++;
	    }
	  }
	}//end of row loop
	if ( nabove_thresh>0 ) {
	  start_nplanes_above_thresh++;
	}
      }//end of plane loop

      if ( crosses_image_boundary==0 ) {
	// start outside image crosses into it
	track_start_boundary = larlitecv::kImageEnd;
	start_crossingname = "ImageEnd";
	start_crosses = true;
      }
      else if ( track_start_dwall<10.0 ) {
	// else has to cross a real boundary. does it?
	start_crosses = true;
      }
      
      // store the track end if it crossed a deetector boundary or the image boundary
      if ( start_crosses ) {
	data.start_pixels.push_back( start_pix );
	data.start_crossingpts.emplace_back( std::move(sce_start) );
	data.start_crossing_flashindex.push_back( flash_index );
	data.start_type.push_back( track_start_boundary );
	data.start_crossing_nplanes_w_charge.push_back( start_nplanes_above_thresh );
	data.tot_true_crossingpoints++;
	data.true_crossingpoints[track_start_boundary]++;
      }
      else {
	track_start_boundary = larlitecv::kUndefined;
	start_crossingname = " [doesnotcross] ";
      }
	

      // ----------------------------------------------------------------------
      // TRACK END
      
      bool end_crosses = false;
      int end_nplanes_above_thresh = 0;
      
      // change end_pix[0] from tick to row
      end_pix[0]   = meta.row( end_pix[0] );

      // check neighborhood around pixel (has charge around?)
      for (int p=0; p<3; p++) {
	int nabove_thresh = 0;
	for (int dr=-endpt_radius; dr<=endpt_radius; dr++) {
	  int srow = end_pix[0]+dr;
	  if ( srow<0 || srow>=(int)meta.rows() )
	    continue;
	  for (int dc=-endpt_radius; dc<=endpt_radius; dc++) {
	    int scol = end_pix[p+1]+dc;
	    if ( scol<0 || scol>=(int)meta.cols() )
	      continue;
	    if ( img_v[p].pixel( srow, scol )>10.0 ) {
	      nabove_thresh++;
	    }
	  }
	}//end of row loop
	if ( nabove_thresh>0 ) {
	  end_nplanes_above_thresh++;
	}
      }//end of plane loop


      if ( crosses_image_boundary==1 ) {
	// end is outside after crossing out of image
	track_end_boundary = larlitecv::kImageEnd;
	end_crossingname = "ImageEnd";
	end_crosses = true;
      }
      else if ( track_end_dwall<10.0 ) {
	// crosses (or gets close enough) to real boundary
	end_crosses = true;
      }

      if ( end_crosses ) {
	data.end_pixels.push_back( end_pix );
	data.end_crossingpts.emplace_back( std::move(sce_end) );
	data.end_crossing_nplanes_w_charge.push_back( end_nplanes_above_thresh );
	data.end_type.push_back( track_end_boundary );
	data.tot_true_crossingpoints++;
	data.true_crossingpoints[track_end_boundary]++;
      }
      else {
	track_end_boundary = larlitecv::kUndefined;
	end_crossingname = " [doesnotcross] ";
      }

      // If we made crossing points, we want to associate the indices of the crossing point info to the mc track
      if ( start_crosses || end_crosses ) {
	data.mctrack_imgendpoint_indices[trackindex].resize(2,-1);
	if ( start_crosses )
	  data.mctrack_imgendpoint_indices[trackindex][0] = data.start_pixels.size()-1;
	if ( end_crosses )
	  data.mctrack_imgendpoint_indices[trackindex][1] = data.end_pixels.size()-1;	  
      }
      
      
      // ----------------------------------------------------------------------------
      // Print summary
      
      if ( printFlashEnds ) {
	std::cout << "  crosses_image_boundary=" << crosses_image_boundary << std::endl;
	std::cout  << "  Found flash: " << start_found_flash << " index=" << flash_index << std::endl;
	
	if ( start_crosses )
	  std::cout << "  Start Boundary Crossing: boundary=" << start_crossingname << " dwall=" << track_start_dwall << " ";
	else
	  std::cout << "  Start did not cross boundary. dwall=" << track_start_dwall << " ";

	std::cout << " row=" << start_pix[0] << " tick=" << meta.pos_y(start_pix[0]) << " (orig=" << orig_start_tick << ")"
		  << " pos=(" << first_step.X() << "," << first_step.Y() << "," << first_step.Z() << ")";
	if ( track_start_boundary==larlitecv::kImageEnd )
	  std::cout << " imgx row=" << data.start_pixels.back()[0]
		    << " imgx pos=(" << data.start_crossingpts.back()[0] << "," << data.start_crossingpts.back()[1] << "," << data.start_crossingpts.back()[2] << ")";
	std::cout << " nplaneswq=" << data.start_crossing_nplanes_w_charge.back()
		  << std::endl;
	
	if ( end_crosses )
	  std::cout << "  End Boundary Crossing: boundary=" << end_crossingname  << " dwall=" << track_end_dwall;
	else
	  std::cout << "  Track End did not cross boundary: dwall=" << track_end_dwall;
	    
	std::cout << " row=" << end_pix[0] << " tick=" << meta.pos_y(end_pix[0]) << " (orig=" << orig_end_tick << ")"
		  << " pos=(" << last_step.X() << "," << last_step.Y() << "," << last_step.Z() << ")";
	if ( track_end_boundary==larlitecv::kImageEnd )
	  std::cout << " imgx row=" << data.end_pixels.back()[0]
		    << " imgx pos=(" << data.end_crossingpts.back()[0] << "," << data.end_crossingpts.back()[1] << "," << data.end_crossingpts.back()[2] << ")";
	std::cout << " nplaneswq=" << data.end_crossing_nplanes_w_charge.back()	<< std::endl;
      }
	
      //if ( start_intime || end_intime ) {
	//intime_cosmics++;
	// if ( start_intime && !start_crosses ) {
	//   std::cout << "start point does not cross boundary: (" << fstart[0] << "," << fstart[1] << "," << fstart[2] << ")"
	// 	    << " true start: (" << track.front().X() << "," << track.front().Y() << "," << track.front().Z() << ") t=" << track.front().T()
	// 	    << " tick = " << getTick( track.front(), ev_trigger->TriggerTime(), &sce )
	// 	    << " dwall=" << track_start_dwall << " intime=" << start_intime
	// 	    << " tick=" << meta.pos_y( start_pix[0] )
	// 	    << " row=" << start_pix[0]
	// 	    << std::endl;
	//   //throw std::runtime_error("start point does not cross boundary?");
	// }
      if ( start_crosses && end_crosses )
	data.true_intime_thrumu++;
      else if ( start_crosses && !end_crosses )
	data.true_intime_stopmu++;
      
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

  int doesTrackCrossImageBoundary( const larlite::mctrack& track, const larcv::ImageMeta& meta, const float trig_time, const larlitecv::SpaceChargeMicroBooNE* psce ) {
    float tick_start = getTick( track.front(), trig_time, psce );
    float tick_end   = getTick( track.back(), trig_time, psce );
    if ( tick_start>meta.min_y() && tick_start<meta.max_y() && tick_end>meta.min_y() && tick_end<meta.max_y() )
      return -1;

    if ( tick_start < meta.min_y() && tick_end > meta.min_y() )
      return 0; // start out -> end in
    else if ( tick_start > meta.min_y() && tick_end < meta.min_y())
      return 1; // start in -> end out
    else if ( tick_start < meta.max_y() && tick_end > meta.max_y())
      return 1; // start in -> end out;
    else if ( tick_start > meta.max_y() && tick_end < meta.max_y() )
      return 0; // start out -> end in

    return -1;
  }

  std::vector<int> getImageBoundaryCrossingPoint( const larlite::mctrack& track, std::vector<float>& crossingpt, const larcv::ImageMeta& meta,
						  const float boundary_tick_buffer, const float trig_time, const larlitecv::SpaceChargeMicroBooNE* psce ) {
    
    if ( doesTrackCrossImageBoundary( track, meta, trig_time, psce )==-1 ) {
      std::stringstream msg;
      msg << __FILE__ << ":" << __LINE__ << " asking for bundary crossing point for track that does not cross the boundary" << std::endl;
      throw std::runtime_error( msg.str() );
    }
    
    crossingpt.resize(3,0);
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    
    float last_tick = getTick( track.front(), trig_time, psce );
    float fpos_last[3] = { (float)track.front().X(), (float)track.front().Y(), (float)track.front().Z() };
    for (int i=1; i<(int)track.size(); i++) {
      const auto& step_now = track[i];
      const auto& step_last = track[i-1];
      float fpos_now[3] = { (float)step_now.X(), (float)step_now.Y(), (float)step_now.Z() };
      float tick_now = getTick( step_now, trig_time, psce );
      float high_tick = tick_now;
      float low_tick  = last_tick;
      if ( low_tick>high_tick ) {
	high_tick = last_tick;
	low_tick = tick_now;
      }
      float boundary_tick = meta.min_y();
      bool crosses_bounds = false;
      if ( low_tick<meta.min_y() && high_tick>meta.min_y() ) {
	crosses_bounds = true;
	boundary_tick = meta.min_y() + boundary_tick_buffer;
      }
      else if ( low_tick<meta.max_y() && high_tick>meta.max_y() ) {
	crosses_bounds = true;
	boundary_tick = meta.max_y() - boundary_tick_buffer;	
      }

      /// go to next step, if not the boundary crossing step
      if ( crosses_bounds ) {

	//std::cout << "found crossing step: tick now=" << tick_now << " last=" << last_tick;
	
	float dir[3] = {0};
	float dirnorm = 0;
	for (int i=0; i<3; i++) {
	  dir[i] = fpos_now[i]-fpos_last[i];
	  dirnorm += dir[i]*dir[i];
	}
	dirnorm = sqrt(dirnorm);
	for (int i=0; i<3; i++) {
	  dir[i] /= dirnorm;
	}
	float dtick = boundary_tick-last_tick;
	float dcm   = cm_per_tick*dtick;
	for (int i=0; i<3; i++) {
	  crossingpt[i] = fpos_last[i] + dcm*dir[i];
	}
	//std::cout << " dtick=" << dtick << " dcm=" << dcm << std::endl;
	// get the image coordinates. Should be within the image now.
	if ( psce ) {
	  std::vector<double> offset = psce->GetPosOffsets( crossingpt[0], crossingpt[1], crossingpt[2] );
	  crossingpt[0] += -offset[0] + 0.7;
	  crossingpt[1] += offset[1];
	  crossingpt[2] += offset[2];
	}
	std::vector<int> crossing_imgcoords = larcv::UBWireTool::getProjectedImagePixel( crossingpt, meta, 3 );
	float finaltick = getTick( step_last, trig_time, psce ) + dtick;
	crossing_imgcoords[0] = meta.row( finaltick );
	return crossing_imgcoords;
      }//end of if crosses

      // if not, update last tick info and move on
      last_tick = tick_now;
      std::memcpy( fpos_last, fpos_now, sizeof(float)*3 );
    }
    // should not get here
    std::stringstream msg;
    msg << __FILE__ << ":" << __LINE__ << " routine should not get here." << std::endl;
    throw std::runtime_error( msg.str() );
    std::vector<int> empty(4,0);
    return empty;
  }
  

  float getTick( const larlite::mcstep& step, const float trig_time, const larlitecv::SpaceChargeMicroBooNE* psce ) {
    std::vector<float> pos(4,0);
    pos[0] = step.T();
    pos[1] = step.X();
    pos[2] = step.Y();
    pos[3] = step.Z();
    return getTick( pos, trig_time, psce );
  }
  

  float getTick( const std::vector<float>& step, const float trig_time, const larlitecv::SpaceChargeMicroBooNE* psce ) {    
    // Function returns the tick time of a MC step point
    // if SCE pointer is null, we do not correct for the space charge
    
    std::vector<double> dpos(3,0);
    if ( psce ) {
      std::vector<double> pos_offset = psce->GetPosOffsets( step[1], step[2], step[3] );
      dpos[0] = step[1] - pos_offset[0] + 0.7;
    }
    else {
      dpos[0] = step[1];
    }
    
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;    
    float tick = ( step[0]*1.0e-3 - (trig_time-4050.0) )/0.5 + dpos[0]/cm_per_tick + 3200.0;
    
    return tick;
  }

  std::vector<float> getFirstStepPosInsideImage( const larlite::mctrack& track, const larcv::ImageMeta& meta, const float trig_time,
						 const bool startAtstart, const float max_step_size, const float fv_border, const larlitecv::SpaceChargeMicroBooNE* psce ) {
    // This function returns the (SCE-corrected) position where a MC track first is inside the image bounds
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;    
    int npts = (int)track.size();

    for ( int ipt=1; ipt<npts; ipt++ ) {

      int thispt = ipt;
      int lastpt = thispt-1;
      if ( !startAtstart ) {
	thispt = npts-1-ipt;
	lastpt = thispt+1;
      }
      
      const auto& this_step = track.at( thispt );
      const auto& last_step = track.at( lastpt );
      
      float dir[3] = { float(this_step.X()-last_step.X()), float(this_step.Y()-last_step.Y()), float(this_step.Z()-last_step.Z()) };
      float dirnorm = 0;
      for (int i=0; i<3; i++) {
	dirnorm += dir[i]*dir[i];
      }
      dirnorm = sqrt(dirnorm);
      if ( dirnorm<1.0e-3 )
	continue;
	
      for (int i=0; i<3; i++)
	dir[i] /= dirnorm;
      int nsteps=dirnorm/max_step_size+1;
      if ( nsteps<= 0 )
	nsteps = 1;
      float stepsize = dirnorm/float(nsteps);
      for (int istep=0; istep<nsteps; istep++) {
	std::vector<float> pos(4,0);
	std::vector<float> pos4v(4,0);
	pos4v[0] = last_step.T();	
	pos[0] = pos4v[1] = last_step.X() + stepsize*float(istep)*dir[0];
	pos[1] = pos4v[2] = last_step.Y() + stepsize*float(istep)*dir[1];
	pos[2] = pos4v[3] = last_step.Z() + stepsize*float(istep)*dir[2];
	
	std::vector<double> offset = psce->GetPosOffsets( pos[0], pos[1], pos[2] );
	std::vector<float> pos_sce(3);
	pos_sce[0] = pos[0]-(float)offset[0]+0.7;
	pos_sce[1] = pos[1]+(float)offset[1];
	pos_sce[2] = pos[2]+(float)offset[2];
	int boundary_type = -1;
	float fdwall = larlitecv::dwall( pos, boundary_type );
	if ( fdwall<fv_border )
	  continue;

	float tick = getTick( pos4v, trig_time, psce );
	if ( tick<meta.min_y()+20.0 || tick>meta.max_y()-20.0 )
	  continue;

	std::cout << " [" << thispt << "/" << npts << ":" << istep << "/" << nsteps << "] tick=" << tick;
	if ( startAtstart )
	  std::cout << " trackstart=(" << track.front().X() << "," << track.front().Y() << "," << track.front().Z() << ")";
	else
	  std::cout << " trackend=(" << track.back().X() << "," << track.back().Y() << "," << track.back().Z() << ")";	  
	std::cout << " truepos=(" << pos[0] << "," << pos[1] << "," << pos[2] << ") ";
	std::cout << " intime pos_sce=(" << pos_sce[0] << "," << pos_sce[1] << "," << pos_sce[2] << ") tick=" << tick << " ";
	pos_sce[0] = (tick-3200.0)*cm_per_tick; // we have to give the apparent-x (relative to the trigger) because we need to test the position in the image	
	std::vector<int> imgcoords;
	try {
	  imgcoords = larcv::UBWireTool::getProjectedImagePixel( pos_sce, meta, 3 );
	  std::cout << " imgcoords=(row=" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")" << std::endl;
	}
	catch (...) {
	  std::cout << std::endl;
	  continue;
	}
	
	return pos;
      }
      
    }//end of track loop

    // didn't find the crossing boundary
    std::vector<float> empty;
    return empty;
  }

  
}
