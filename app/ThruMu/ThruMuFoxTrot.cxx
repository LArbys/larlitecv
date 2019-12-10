#include "ThruMuFoxTrot.h"

#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  // --------------------------------------------------------------------------
  // ThruMu Lead. Choose position most consistent with straight-line after
  // trying to remove space-charge effect.
  bool ThruMuFoxTrotLead::_chooseBestSegment_( const FoxStep& current, std::vector<Segment3D_t>& seg3d_v, const FoxTrack& past,
					       const FoxTrotTrackerAlgoConfig& config, const std::vector<larcv::Image2D>& stepped_v, int& best_seg_idx ) {

    float closest_dist2realline = 1.0e9;
    best_seg_idx = -1;
    for ( int iseg=0; iseg<(int)seg3d_v.size(); iseg++ ) {
      Segment3D_t& seg = seg3d_v[iseg];

      // don't go backwards
      float segdist = 0.;
      float dirnorm = 0.;
      float segdir[3] = {0};
      for (int i=0; i<3; i++) {
        segdir[i] = seg.end[i]-seg.start[i];
        segdist += (seg.end[i]-seg.start[i])*(seg.end[i]-seg.start[i]);
        dirnorm += current.dir()[i]*current.dir()[i];
      }
      segdist = sqrt(segdist);
      dirnorm = sqrt(dirnorm);

      float cosdir = 0.;
      for ( int i=0; i<3; i++ ){
        cosdir += current.dir()[i]*segdir[i]/(segdist*dirnorm);
      }


      // get original pos
      std::vector<double> segend(3,0);
      for (int i=0; i<3; i++)
        segend[i] = seg.end[i];
      segend[0] += m_shiftx;
      std::vector<double> doriginal_pos = m_sce.getOriginalPos( segend );

      double dist2realline = m_geoalgo.SqDist( real_line3d, doriginal_pos );

      std::cout << "seg3d #" << iseg << " dist2realline = " << dist2realline << " cosdir=" << cosdir << std::endl;

      if ( cosdir>0.0 && dist2realline < closest_dist2realline ) {
        closest_dist2realline = dist2realline;
        best_seg_idx = iseg;
      }

    }//end of segment loop

    if ( best_seg_idx>=0 )
      return true;

    return false;
  }

  double ThruMuFoxTrotLead::getDist2RealLine( const std::vector<double>& pt ) {
    return sqrt( m_geoalgo.SqDist( real_line3d, pt ) );
  }

  std::vector<double> ThruMuFoxTrotLead::getClosestPointOnRealLine( const std::vector<double>& pt ) {
    return m_geoalgo.ClosestPt( pt, real_line3d );
  }

  const std::vector<double>& ThruMuFoxTrotLead::getRealLineDir() const {
    return m_realdir;
  }

  std::vector<double> ThruMuFoxTrotLead::getSCEAppliedPos( const std::vector<double>& pos ) {
    std::vector<double> offset = m_sce_forward.GetPosOffsets( pos[0], pos[1], pos[2] );
    std::vector<double> scepos(3,0);
    scepos[0] = pos[0] - offset[0] + 0.17;
    scepos[1] = pos[1] + offset[1];
    scepos[2] = pos[2] + offset[2];
    return scepos;
  }

  void ThruMuFoxTrotLead::setEndPoints( const std::vector<float>& start, const std::vector<float>& end, const float shiftx ) {

    m_shiftx = shiftx;

    std::vector<float> shift_start = start;
    shift_start[0] += m_shiftx;
    std::vector<float> shift_end   = end;
    shift_end[0]   += m_shiftx;

    image_line3d.Pt1( shift_start[0], shift_start[1], shift_start[2] );
    image_line3d.Pt2( shift_end[0],   shift_end[1],   shift_end[2] );

    std::vector<float> original_start = m_sce.getOriginalPos( shift_start );
    std::vector<float> original_end   = m_sce.getOriginalPos( shift_end );
    real_line3d.Pt1( original_start[0], original_start[1], original_start[2] );
    real_line3d.Pt2( original_end[0], original_end[1], original_end[2] );

    double norm = 0.;
    m_realdir.resize(3,0);
    for (int i=0; i<3; i++) {
      double dx = original_end[i] - original_start[i];
      std::cout << "( original_end[" << i << "]= " << original_end[i] << " - original_start[" << i << "]=" << original_start[i] << ") =" << dx << std::endl;
      m_realdir[i] = dx;
      norm += dx*dx;
    }
    norm = sqrt(norm);
    for (int i=0; i<3; i++)
      m_realdir[i] /= norm;
    std::cout << "set thrumufoxtrot real line: dir=(" << m_realdir[0] << "," << m_realdir[1] << "," << m_realdir[2] << ")" << std::endl;
  }

  void ThruMuFoxTrotLead::setEndPoints( const std::vector<double>& start, const std::vector<double>& end, const double shiftx ) {
    std::vector<float> fstart(3,0.0);
    std::vector<float> fend(3,0.0);
    for (int i=0; i<3; i++) {
      fstart[i] = start[i];
      fend[i] = end[i];
    }
    setEndPoints( fstart, fend, shiftx );
  }

  // --------------------------------------------------------------------------
  // THRUMU FOXTROT

  ThruMuFoxTrot::ThruMuFoxTrot( const ThruMuFoxTrotConfig& cfg )
    : m_config(cfg), m_tracker(cfg.foxtrotalgo_cfg) {
      if (m_config.use_thrumu_lead )
        m_tracker.setUserLead( &m_lead );
  }


  BMTrackCluster3D ThruMuFoxTrot::findThruMuTrack( const BoundarySpacePoint& pt_a, const BoundarySpacePoint& pt_b,
						   const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v ) {
    BMTrackCluster3D bmtrack;
    return bmtrack;
  }

  FoxTrack ThruMuFoxTrot::findThruMuFoxTrack( const BoundarySpacePoint& pt_a, const BoundarySpacePoint& pt_b,
					      const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v ) {
    // we give the end points to the FoxTrot lead

    std::cout << "=========== THRUMU FOX TROT =======================" << std::endl;

    // calculate shift
    float shiftx = 0.;
    if ( pt_a.type()==larlitecv::kAnode ) {
      shiftx = -pt_a.pos()[0]; // should be zero
    }
    else if ( pt_b.type()==larlitecv::kAnode ) {
      shiftx = -pt_b.pos()[0]; // should be zero
    }
    else if ( pt_a.type()==larlitecv::kCathode ) {
      shiftx = 258.0 - pt_a.pos()[0];
    }
    else if ( pt_b.type()==larlitecv::kCathode ) {
      shiftx = 258.0 - pt_b.pos()[0];
    }
    std::cout << "Shift-x: " << shiftx << std::endl;

    // do we flip the start and end point? we favor the use of either anode or cathode, because those tracks have a known x position by definition.
    const BoundarySpacePoint* start = &pt_a;
    const BoundarySpacePoint* end   = &pt_b;
    if ( (pt_b.type()==larlitecv::kCathode || pt_b.type()==larlitecv::kAnode ) && ( pt_a.type()!=larlitecv::kCathode && pt_a.type()!=larlitecv::kAnode ) ) {
      start = &pt_b;
      end   = &pt_a;
    }

    std::cout << "Starting point: (" << (*start).pos()[0] << "," << (*start).pos()[1] << "," << (*start).pos()[2] << ")" << std::endl;
    std::cout << "Ending point: (" << (*end).pos()[0] << "," << (*end).pos()[1] << "," << (*end).pos()[2] << ")" << std::endl;

    m_lead.setEndPoints( (*start).pos(), (*end).pos(), shiftx );

    FoxTrack track;
    m_step_types.clear();
    bool finished = false;

    BoundarySpacePoint sp_start = (*start);
    double last_dist2goal = 1.0e9;
    int num_steps_away_from_goal = 0;
    int num_iters = 0;

    while ( !finished ) {

      std::cout << "NUM ITERS: " << num_iters << std::endl;

      // we use the fox trot tracker algo to follow lines of charge from start to end
      FoxTrack trot = m_tracker.followTrack( img_v, badch_v, badch_v, sp_start );

      // we get to the end?
      // we probably have to scan backwards
      double dist2end = distanceToTheEnd( trot.back(), (*end) );
      std::cout << "Dist2End=" << dist2end << std::endl;
      if ( dist2end<m_config.endpoint_radius ) {
        std::cout << "holy shit, got to the end point. dist2end=" << dist2end << std::endl;
        // fill the track
        for ( auto& step : trot ) {
          track.push_back( step );
          m_step_types.push_back(0);
        }
        break;
      }


      cleanTrackEnd( trot, 0.8 );

      // fill the track
      for ( auto& step : trot ) {
        track.push_back( step );
        m_step_types.push_back(0);
      }



      // split the track
      std::vector< T3DCluster > t3dtracks = makeRealAndImagePaths( trot );
      std::vector< std::vector<T3DCluster> > t3d_v = breakPathsIntoStraightSegments( t3dtracks );

      // we see if we are going to append to the latest split
      if ( m_track_splits.size()==0 ) {
        // nothing to merge against
        m_track_splits.emplace_back( std::move(t3d_v.front()) );
        std::cout << "PUT IN THE FIRST SPLIT" << std::endl;
      }
      else {
        // test if we should merge
        T3DCluster& past_realtrack    = m_track_splits.back()[1];
        T3DCluster& current_realtrack = t3d_v.front()[1];
        std::vector<int> whichends(2,0);
        whichends[0] = 0;
        whichends[1] = 1;
        double closest_dist = 0;
        bool letsmerge = m_pcmerge.shouldWeEndPointMerge( past_realtrack, current_realtrack, closest_dist, whichends );
        if ( letsmerge ) {
          past_realtrack.append( current_realtrack );
          T3DCluster& past_imagetrack = m_track_splits.back()[0];
          T3DCluster& past_realtrack  = t3d_v.front()[0];
          past_imagetrack.append( past_realtrack );
          std::cout << "MERGE THE CURRENT SPLIT" << std::endl;
        }
        else {
          m_track_splits.emplace_back( std::move(t3d_v.front()) );
          std::cout << "UPDATE THE CURRENT SPLIT" << std::endl;
        }
      }

      if ( t3d_v.size()==2 ) {
        // regardless if we merged the split, we need to put the latter split
        m_track_splits.emplace_back( std::move(t3d_v.back()) );
        std::cout << "ADD THE LATEST SPLIT" << std::endl;
      }
      std::cout << "After FoxTrotTracker, number of splits=" << m_track_splits.size() << std::endl;
      std::cout << "track point after trot: (" << trot.back().pos()[0] << "," << trot.back().pos()[1] << "," << trot.back().pos()[2] << ") dist2end=" << dist2end << " cm" << std::endl;

      // update the FoxTrotLead
      std::vector< T3DCluster >& cluster_splits = m_track_splits.back();
      std::vector<float> current_start(3);
      for (int i=0; i<3; i++)
        current_start[i] = cluster_splits[0].getPath().back()[i];
      m_lead.setEndPoints( current_start, (*end).pos(), shiftx ); // we set the image-view end points

      // Hypothesis track
      FoxTrack hypotrack;
      FoxStep* current_hypo = &(trot.back());
      for (int istep=0; istep<20; istep++) {
        FoxStep step = makeHypothesisStep( *current_hypo, 3.0 );
        std::cout << "thrumu hypothesis foxstep=(" << step.pos()[0] << "," << step.pos()[1] << "," << step.pos()[2] << ")" << std::endl;
        hypotrack.push_back( step );
        current_hypo = &(hypotrack.back());
        //track.emplace_back( std::move(step) );
        if ( doesHypothesisStepSeeCharge( step, img_v, badch_v ) ) {
          sp_start.setZY( step.pos()[2], step.pos()[1] );
          sp_start.setX( step.pos()[0] );
          break;
        }
      }

      std::cout << "Next starting point: (" << sp_start.pos()[0] << "," << sp_start.pos()[1] << "," << sp_start.pos()[2] << ")" << std::endl;

      std::vector< std::vector<double> > image_hypo;
      std::vector< std::vector<double> > real_hypo;
      for ( auto& step : hypotrack ) {
        track.push_back( step );
        m_step_types.push_back(1);

        std::vector<double> imagepos(3);
        std::vector<double> realpos(3);
        for (int i=0; i<3; i++) {
          imagepos[i] = step.pos()[i];
          realpos[i]  = step.pos()[i];
        }

        image_hypo.push_back( imagepos );
        real_hypo.push_back( m_lead.getSCECorrectedPos( realpos ) );

      }

      if ( image_hypo.size()>0 ) {
        // append to the current T3D split
        T3DCluster t3d_image_hypo( image_hypo );
        T3DCluster t3d_real_hypo( real_hypo );
        m_track_splits.back()[0].append( t3d_image_hypo );
        m_track_splits.back()[1].append( t3d_real_hypo );
      }

      //break;

      double dist2goal = distanceToTheEnd( track.back(), (*end) );
      if ( last_dist2goal < dist2goal )
        num_steps_away_from_goal++;
      std::cout << "distance to goal=" << dist2goal << ".  we've taken " << num_steps_away_from_goal << " iterations away from our goal. "
                << " iteration=" << num_iters << " of " << m_config.maxiters << std::endl;
      last_dist2goal = dist2goal;

      if ( dist2goal<m_config.endpoint_radius || num_steps_away_from_goal>=2 || num_iters>=m_config.maxiters)
      //if ( dist2goal<m_config.endpoint_radius || num_steps_away_from_goal>=2 || num_iters>=2 )
        finished = true;

      std::cout << "END OF ITERATION #" << num_iters << std::endl;
      std::cin.get();
      num_iters++;
      break;
      //*/
    }

    return track;
  }

  FoxStep ThruMuFoxTrot::makeHypothesisStep( const FoxStep& current_step, const double step_size ) {
    // we make a hypothesis point using space charge corrections and the straight-line thrumu hypothesis
    std::vector< double > dpos(3,0);
    for (int i=0; i<3; i++)
      dpos[i] = current_step.pos()[i];
    std::vector<double> sce_corrected_pt = m_lead.getSCECorrectedPos( dpos );
    double dist2realline = m_lead.getDist2RealLine( sce_corrected_pt );
    std::vector<double> pt_on_realline = m_lead.getClosestPointOnRealLine( sce_corrected_pt );

    std::cout << "pt_on_realline=(" << pt_on_realline[0] << "," << pt_on_realline[1] << "," << pt_on_realline[2] << ")" << std::endl;

    const std::vector<double>& realdir = m_lead.getRealLineDir();
    std::cout << "realline dir=(" << realdir[0] << "," << realdir[1] << "," << realdir[2] << ")" << std::endl;

    std::vector<double> realpt2scept_dir(3,0);
    double norm = 0;
    for ( int i=0; i<3; i++) {
      realpt2scept_dir[i] = sce_corrected_pt[i] - pt_on_realline[i];
      norm += realpt2scept_dir[i]*realpt2scept_dir[i];
    }
    norm = sqrt(norm);
    for (int i=0; i<3; i++) {
      realpt2scept_dir[i] /= norm;
    }

    if ( dist2realline<1.0 ) // ignore the noise
      dist2realline = 0.0;

    std::vector<double> next_real_step(3,0);
    for (int i=0; i<3; i++) {
      next_real_step[i] = pt_on_realline[i] + realdir[i]*step_size + realpt2scept_dir[i]*dist2realline;
    }

    // space charge transform
    std::vector<double> next_sce_pt = m_lead.getSCEAppliedPos( next_real_step );

    std::vector<double> next_dir(3);
    double next_norm = 0.;
    for (int i=0; i<3; i++) {
      next_dir[i] = next_sce_pt[i] - current_step.pos()[i];
      next_norm += next_dir[i]*next_dir[i];
    }
    next_norm = sqrt(next_norm);
    for (int i=0; i<3; i++)
      next_dir[i] /= next_norm;

    FoxStep nextfoxstep( next_sce_pt, next_dir );
    return nextfoxstep;
  }

  bool ThruMuFoxTrot::doesHypothesisStepSeeCharge( FoxStep& hypostep, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v ) {
    std::vector<int> image_coords = larcv::UBWireTool::getProjectedImagePixel( hypostep.pos(), img_v.front().meta(), img_v.size() );
    const larcv::ImageMeta& meta = img_v.front().meta();
    if ( image_coords[0]<0 || image_coords[0]>meta.rows() ) {
      return false;
    }
    for (int ic=1; ic<(int)image_coords.size(); ic++) {
      if ( image_coords[ic]<0 || image_coords[ic]>meta.cols() )
	return false;
    }

    int numplanes_w_charge = 0;
    int numplanes_chargeorbadch = 0;
    int row = image_coords[0];
    for (size_t p=0; p<img_v.size(); p++) {
      int col = image_coords[p+1];
      bool hascharge = false;
      bool hasbadch  = false;
      for (int dr=-m_config.hit_neighborhood; dr<=m_config.hit_neighborhood; dr++) {
	int r = row+dr;
	if ( r<0 || r>(int)meta.rows() )
	  continue;
	for (int dc=-m_config.hit_neighborhood; dc<=m_config.hit_neighborhood; dc++) {
	  int c = col+dc;
	  if ( c<0 || c>(int)meta.cols() )
	    continue;

	  if ( img_v[p].pixel( r, c )>m_config.pixel_thresholds[p] )
	    hascharge = true;

	  if ( badch_v[p].pixel( r, c )>0 )
	    hasbadch = true;
	} //end of dc loop
      }//end of dr loop

      if ( hascharge )
	numplanes_w_charge++;
      if ( hascharge || hasbadch )
	numplanes_chargeorbadch++;
    }//end of plane loop

    if ( numplanes_w_charge>=(int)img_v.size()-1 && numplanes_chargeorbadch>=(int)img_v.size() )
      return true;

    return false;
  }

  double ThruMuFoxTrot::distanceToTheEnd( FoxStep& current, const BoundarySpacePoint& endpoint ) {
    double dist = 0.;
    for (int i=0; i<3; i++) {
      dist += ( current.pos()[i]-endpoint.pos()[i] )*( current.pos()[i]-endpoint.pos()[i] );
    }
    dist = sqrt(dist);

    return dist;
  }

  void ThruMuFoxTrot::cleanTrackEnd( FoxTrack& track, float cos_requirement ) {
    int steps_before = track.size();
    int end = (int)track.size()-5;
    if ( end<0 )
      end = 0;

    for (int istep=(int)track.size()-1; istep>=(end+1); istep-- ) {
      float stepcos = 0.;
      for (int i=0; i<3; i++) {
	stepcos += track[istep].dir()[i]*track[istep-1].dir()[i];
      }
      std::cout << "step #" << istep << ": " << stepcos << std::endl;
    }

    int remove_at = -1;
    for (int istep=(int)track.size()-1; istep>=1; istep-- ) {
      float stepcos = 0.;
      for (int i=0; i<3; i++) {
	stepcos += track[istep].dir()[i]*track[istep-1].dir()[i];
      }
      std::cout << "step #" << istep << ": " << stepcos << std::endl;
      if ( stepcos<cos_requirement ) {
	remove_at = istep;
      }
      else if ( steps_before-istep>=5 )
	break;
    }

    if ( remove_at>=0 ) {
      for ( int istep=(int)track.size()-1; istep>=remove_at; istep-- )
	track.pop_back();
    }

    std::cout << "Steps removed: " << steps_before - (int)track.size() << std::endl;
  }

  std::vector< T3DCluster > ThruMuFoxTrot::makeRealAndImagePaths( FoxTrack& trot ) {

    std::vector< std::vector<double> > image_path;
    std::vector< std::vector<double> > real_path;
    for ( auto& step : trot ) {

      // apparent path (in the image)
      std::vector<double> pos(3,0);
      for (int i=0; i<3; i++)
	pos[i] = step.pos()[i];
      image_path.push_back( pos );

      // real path
      std::vector<double> real_pos = m_lead.getSCECorrectedPos( pos );
      real_path.push_back( real_pos );
    }

    T3DCluster t3d_image( image_path );
    T3DCluster t3d_real( real_path );

    std::vector< T3DCluster > vec;
    vec.emplace_back( std::move(t3d_image) );
    vec.emplace_back( std::move(t3d_real) );

    return vec;
  }

  std::vector< std::vector< T3DCluster > > ThruMuFoxTrot::breakPathsIntoStraightSegments( std::vector<T3DCluster>& paths ) {

    // this code looks to split track based on defect
    T3DCluster& image_path = paths[0];
    T3DCluster& real_path  = paths[1];

    ::larlite::geoalgo::Line start2end( real_path.getPath().front(), real_path.getPath().back() );

    double maxdist = 0;
    int max_istep = 0;
    for (int istep=0; istep<(int)real_path.getPath().size(); istep++) {
      auto& realpos = real_path.getPath()[istep];
      double dist = m_geoalgo.SqDist( start2end, realpos );
      dist = sqrt(dist);
      if ( dist>maxdist ){
	maxdist = dist;
	max_istep = istep;
      }
    }

    // split this
    std::cout << "thrumufoxtrot defect dist=" << maxdist << " @ "
	      << "(" << real_path.getPath()[max_istep][0] << "," << real_path.getPath()[max_istep][1] << "," << real_path.getPath()[max_istep][2] << ")"
	      << std::endl;

    std::vector< std::vector< T3DCluster > > splits;
    if ( maxdist>10.0 ) {

      std::vector< std::vector<double> > imagepath1;
      std::vector< std::vector<double> > realpath1;
      for (int istep=0; istep<=max_istep; istep++) {
	imagepath1.push_back( image_path.getPath()[istep] );
	realpath1.push_back( real_path.getPath()[istep] );
      }
      std::vector< T3DCluster > split1;
      T3DCluster imagetrack1( imagepath1 );
      T3DCluster realtrack1( realpath1 );
      split1.emplace_back( std::move(imagetrack1) );
      split1.emplace_back( std::move(realtrack1) );

      std::vector< std::vector<double> > imagepath2;
      std::vector< std::vector<double> > realpath2;
      for (int istep=max_istep+1; istep<(int)real_path.getPath().size(); istep++) {
	imagepath2.push_back( image_path.getPath()[istep] );
	realpath2.push_back( real_path.getPath()[istep] );
      }
      std::vector< T3DCluster > split2;
      T3DCluster imagetrack2( imagepath2 );
      T3DCluster realtrack2( realpath2 );
      split2.emplace_back( std::move(imagetrack2) );
      split2.emplace_back( std::move(realtrack2) );


      splits.emplace_back( std::move(split1) );
      splits.emplace_back( std::move(split2) );
    }
    else {
      // pass back
      splits.emplace_back( std::move(paths) );
    }

    return splits;
  }

}
