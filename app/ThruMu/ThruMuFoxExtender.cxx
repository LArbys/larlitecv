#include "ThruMuFoxExtender.h"

#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  bool ThruMuFoxExtender::extendTrack( std::vector<std::vector<double> >& track, const std::vector<larcv::Image2D>& img_v,
				       const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v ) {

    std::vector<float> startpoint(3,0);
    std::vector<float> startdir(3,0);
    std::vector<float> endpoint(3,0);
    std::vector<float> enddir(3,0);
    int npts = track.size();
    float norm = 0.;
    float norm_start = 0.;
    for (int i=0; i<3; i++) {
      endpoint[i] = track[npts-1][i];
      enddir[i] = endpoint[i]-track[npts-2][i];
      norm += enddir[i]*enddir[i];
      startpoint[i] = track[0][i];
      startdir[i]   = track[0][i]-track[1][i];
      norm_start += startdir[i]*startdir[i];
    }
    norm = sqrt(norm);
    norm_start = sqrt(norm_start);
    for (int i=0; i<3; i++) {
      enddir[i] /= norm;
      startdir[i] /= norm_start;
    }

    FoxTrack track_back  = extendFromPoint( endpoint, enddir, img_v, badch_v, tagged_v );
    FoxTrack track_front = extendFromPoint( startpoint, startdir, img_v, badch_v, tagged_v );

    std::vector<float> thresholds(3,10.0);
    int hit_neighborhood = 2;
    bool extended_back  = appendExtension( kBack, track_back, track, 0.3, img_v, badch_v, tagged_v, thresholds, hit_neighborhood );
    bool extended_front = appendExtension( kFront, track_front, track, 0.3, img_v, badch_v, tagged_v, thresholds, hit_neighborhood );

    return extended_front | extended_back;

  }

  FoxTrack ThruMuFoxExtender::extendFromPoint( const std::vector<float>& pos, const std::vector<float>& dir, const std::vector<larcv::Image2D>& img_v,
					       const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v ) {

    FoxTrotTrackerAlgo algo( m_config );

    // make a fake boundary point
    BoundarySpacePoint sp( larlitecv::kUndefined, pos, dir, img_v.front().meta() );
    FoxTrack track = algo.followTrack( img_v, badch_v, tagged_v, sp, dir );

    return track;
  }

  bool ThruMuFoxExtender::appendExtension( ThruMuFoxExtender::ExtEnd_t end, FoxTrack& extension,
					   std::vector<std::vector<double> >& track, const float max_step_size, const std::vector<larcv::Image2D>& img_v,
					   const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
					   const std::vector<float>& thresholds, const int hit_neighborhood ) {
    // evaluate the quality of the extension: this should probably be the extender's function

    // if the extension is less than 2 points, don't append anything
    int nsteps = extension.size();
    if (nsteps<2) {
      return false;
    }

    // we make a list of extension points that go in the right direction (away from the track)
    std::vector<double> startpoint(3,0);
    std::vector<double> endpoint(3,0);
    for (int i=0; i<3; i++) {
      endpoint[i] = track.back()[i];
      startpoint[i] = track.front()[i];
    }
    float lastdist2goal = 0.;
    float lastdist2start = 0.;
    std::vector< std::vector<double> > extension_points;
    extension_points.reserve( nsteps+1 );

    // first point of the extension are the ends of the track
    if ( end==kBack ) {
      extension_points.push_back( track.back() );
    }
    else {
      extension_points.push_back( track.front() );
    }

    // loop over the extension
    for ( int istep=1; istep<nsteps; istep++ ) {
      float dist2goal = 0.;
      float dist2start = 0.;
      for (int i=0; i<3; i++) {
        float dgoal  =  extension[istep].pos()[i]-track.back()[i];
        float dstart =  extension[istep].pos()[i]-track.front()[i];
        dist2goal  += dgoal*dgoal;
        dist2start += dstart*dstart;
      }
      dist2goal  = sqrt(dist2goal);
      dist2start = sqrt(dist2start);
      // probably should be heading to the walls?
      if ( lastdist2goal<dist2goal || lastdist2start<dist2start ) {
        // extend
        std::vector<double> dextpos(3,0);
        for (int i=0; i<3; i++)
	  dextpos[i] = extension[istep].pos()[i];
        extension_points.push_back( dextpos );
        lastdist2start = dist2start;
        lastdist2goal  = dist2goal;
      }
      else {
        break;
      }
    }// loop over extension track

    // if we didn't add anything, return
    if (extension_points.size()<=1)
      return false;

    // now add points taking a max step of max_step_size
    std::vector<std::vector<double> > pts_for_track;

    // loop through extension points determined
    int next = extension_points.size();

    int step_wo_appending = 0; // we stop the extension if we have a certain number of points with charge
    const int stop_after_nnonappend_step = 3;

    for (int iext=1; iext<next; iext+=1) {
      double dir[3] = {0};
      std::vector<double>& pos1 = extension_points[iext-1];
      std::vector<double>& pos2 = extension_points[iext];
      double norm = 0.;
      for (int i=0; i<3; i++) {
	dir[i] = pos2[i]-pos1[i];
	norm += dir[i]*dir[i];
      }
      norm = sqrt(norm);
      for (int i=0; i<3; i++)
	dir[i] /= norm;

      // fox trot steps large, so we break it into substeps
      int nsteps = norm/max_step_size;
      float stepsize = norm/nsteps;
      for (int istep=0; istep<nsteps; istep++) {
	std::vector<float> steppos(3,0);
	std::vector<double> dsteppos(3,0);
	for (int i=0; i<3; i++) {
	  steppos[i]  = pos1[i] + dir[i]*istep*stepsize;
	  dsteppos[i] = pos1[i] + dir[i]*istep*stepsize;
	}
	// project into image and determine if all 3 planes saw charge (or dead channel)
	std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( steppos, img_v.front().meta(), img_v.size() );
	int planes_w_hit = 0;
	for (size_t p=0; p<img_v.size(); p++) {
	  bool hashit = false;
	  for (int dr=-hit_neighborhood; dr<=hit_neighborhood; dr++) {
	    int row = imgcoords[0]+dr;
	    if ( row<0 || row>=(int)img_v[p].meta().rows() ) continue;
	    for (int dc=-hit_neighborhood; dc<=hit_neighborhood; dc++) {
	      int col = imgcoords[p+1]+dc;
	      if ( col<0 || col>=(int)img_v[p].meta().cols() ) continue;
	      if ( img_v[p].pixel(row,col)>thresholds[p] || badch_v[p].pixel(row,col)>0 )
		hashit = true;
	      if (hashit)
		break;
	    }
	    if (hashit)
	      break;
	  }
	  if ( hashit )
	    planes_w_hit++;
	}//end of plane loop

	// determine if we append
	if ( planes_w_hit>=3 ) {
	  // extend
	  pts_for_track.push_back(dsteppos);
	  step_wo_appending = 0; // reset counter
	}
	else {
	  step_wo_appending++;
	}

	// stop after
	if (step_wo_appending>=stop_after_nnonappend_step)
	  break;
      }//end of subsetp loop
      if (step_wo_appending>=stop_after_nnonappend_step)//ed
	break;//ed
    }//end of extension point loop

    // add points to the track
    if ( end==kBack ) {
      // adding to the back is easy
      for (auto& pt : pts_for_track ) {
	track.push_back( pt );
      }
    }
    else {
      // extend the front: harder. loop backwards through extension, then append original track. replace original
      std::vector<std::vector<double> > newtrack;
      newtrack.reserve( track.size()+pts_for_track.size() );
      for (int i=pts_for_track.size()-1; i>=0; i--) {
	newtrack.push_back( pts_for_track[i] );
      }
      for (int i=0; i<track.size(); i++) {
	newtrack.push_back( track[i] );
      }
      track.clear();
      track.reserve( newtrack.size() );
      for (int i=0; i<newtrack.size(); i++)
	track.push_back( newtrack[i] );
    }//end of if front

    return true;
  }
}
