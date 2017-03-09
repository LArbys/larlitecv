#include "StopMuTracker.h"

#include <stdio.h>
#include <time.h>
#include <sstream>

#ifndef __CINT__
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif
#endif

// ROOT
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

// larcv
#include "UBWireTool/UBWireTool.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

#include "StopMuSkeleton.h"

namespace larlitecv {

  StopMuTracker::StopMuTracker( const StopMuTrackerConfig& cfg, const std::vector<larcv::Image2D>& img_v,
				const std::vector<larcv::Image2D>& thrumu_v,
				const std::vector< std::vector< const larcv::Pixel2D* > >& candidate_stopmu_points,
				int verbosity ) : m_config(cfg) {

    setVerbosity(m_config.verbosity);
				
    
    // skeletonize image
    StopMuSkeleton skeleton_op;

    size_t nplanes = m_config.nplanes;

    for (int p=0; p<nplanes; p++) {
      larcv::Image2D skel = skeleton_op.skeletonize( img_v.at(p), m_config.pixel_thresholds[p], m_config.skeleton_kernel_size );
      skel_v.emplace_back( std::move(skel) );
    }

    // with candidate points provided, we make an image marking neighboring pixels.
    // this is to protect these pixels from thrumu-pixel masking step which follows
    std::vector<larcv::Image2D> candidate_pixels;
    for (size_t p=0; p<nplanes; p++) {
      larcv::Image2D img( img_v.at(p).meta() );
      img.paint(0.0);
      for ( auto& pix_v : candidate_stopmu_points ) {
	int row=pix_v.at(p)->Y();
	int col=pix_v.at(p)->X();

	for ( int r=-m_config.start_point_pixel_neighborhood; r<=m_config.start_point_pixel_neighborhood; r++) {
	  int rtest = row+r;
	  if ( rtest<0 || rtest>=(int)img_v.at(p).meta().rows() ) continue;
	  for (int c=-m_config.start_point_pixel_neighborhood; c<=m_config.start_point_pixel_neighborhood; c++) {
	    int ctest = col+c;
	    if ( ctest<0 || ctest>=(int)img_v.at(p).meta().cols() ) continue;
	    if ( img_v.at(p).pixel( rtest, ctest )>m_config.pixel_thresholds[p] )
	      img.set_pixel( rtest, ctest, 1 );
	  }
	}
      }//end of loop over candidate points
      candidate_pixels.emplace_back( std::move(img) );
    }//end of loop over planes
    
    // cluster skeleton pixel, but mask thru-mu pixels
    time_t cluster_start = time(NULL);
    
    for (size_t p=0; p<3; p++) {
      dbscan::dbPoints data;
      const larcv::Image2D& skelimg = skel_v.at(p);
      const larcv::Image2D& candpix = candidate_pixels.at(p);
      for (size_t r=0; r<skelimg.meta().rows(); r++) {
	for (size_t c=0; c<skelimg.meta().cols(); c++) {
	  if ( skelimg.pixel(r,c)==0 ) continue; // not skeleton
	  if ( thrumu_v.at(p).pixel(r,c)>0 && candpix.pixel(r,c)==0 ) continue; // mask thrumu if not near candidate pixel
	  std::vector<double> point(2);
	  point[0] = c; // X 
	  point[1] = r; // Y
	  data.emplace_back( point );
	}
      }
      dbscan::DBSCANAlgo dbalgo;
      dbscan::dbscanOutput cluster = dbalgo.scan( data, m_config.dbscan_cluster_minpoints, m_config.dbscan_cluster_radius );
      m_imghits.emplace_back( data );
      m_clusters.emplace_back( cluster );
      //std::cout << "number of clusters on plane " << p << ": " << m_clusters.at(p).clusters.size() << std::endl;
    }//loop over clusters
    
    time_t cluster_finished = time(NULL);
    double dt_clusters = cluster_finished-cluster_start;
    if ( verbosity>1 )
      std::cout << "clustered in " << dt_clusters << " seconds." << std::endl;
    
    // initialize hit list
    for (size_t p=0; p<m_config.nplanes; p++) {
      current_hit[p] = -1;
    }
#ifdef USE_OPENCV
    // for debug
    // for (int p=0; p<3; p++) {
    //   cv::Mat imgmat = larcv::as_mat_greyscale2bgr( skel_v.at(p), 0, 1.0);
    //   std::stringstream ss;
    //   ss << "skel_p" << p << ".jpg";
    //   cv::imwrite( ss.str().c_str(), imgmat );
    // }
#endif
  }

  // ------------------------------------------------------------------------------------------------
  // *** Tracker Loop Functions ****
  // ------------------------------------------------------------------------------------------------

  void StopMuTracker::makeProposedPos( const std::vector<float>& currentpos, const std::vector<float>& currentdir, std::vector<float>& proposedpos, const float stepsize ) {
    if ( currentpos.size()!=3 && currentdir.size()!=3 ) {
      throw std::runtime_error("StopMuTracker::makeProposedPos length of input pos and dir vectors not 3");
    }
    float dirnorm = 0.;
    proposedpos.resize(3,0.0);
    for (int i=0; i<3; i++) {
      dirnorm += currentdir[i]*currentdir[i];
      proposedpos[i] = currentpos[i] + currentdir[i]*stepsize;
    }
    dirnorm = sqrt(dirnorm);
    if ( std::fabs(dirnorm-1.0)>1.0e-4 ) {
      throw std::runtime_error("StopMuTracker::makeProposedPos direction vector was not normalized to 1");
    }
  }

  void StopMuTracker::imagePositions( const std::vector<float>& currentpos, int& tick, std::vector<int>& wid ) {
    // there is a problem here. this is a static function. so i can't set the trigger time below.  why is this static?
    wid.resize(3,-1); // wire position
    double dpos[3];
    for (int p=0; p<3; p++) dpos[p] = currentpos[p];
    for (int p=0; p<3; p++) {
      wid[p] = (int)::larutil::Geometry::GetME()->WireCoordinate( dpos, p );
      if ( wid[p]<0 ) wid[p] = 0;
      if ( wid[p]>(int)larutil::Geometry::GetME()->Nwires(p) ) wid[p] = (int)larutil::Geometry::GetME()->Nwires(p)-1;
    }
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*(0.5); // [cm/usec]*[usec/tick]
    tick = (int)(currentpos[0]/cm_per_tick) + 3200;
  }

  void StopMuTracker::getClosestHitsInPlane( const int clusterid, const std::vector<int>& test_pos, 
					     const ::dbscan::dbPoints& src_data, const ::dbscan::dbscanOutput& cluster_info, const larcv::ImageMeta& meta, 
					     std::vector< std::pair<int,double> >& hitlist ) {
    // a wrapper around dbscanOutput::closestHitsInCluster
    // have to convert the 2D pixel col/row information into double test_point
    std::vector<double> dtest_pos(test_pos.size(),0);
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*m_config.usec_per_tick;
    float cm_per_wire = m_config.cm_per_wire;
    
    dtest_pos[0] = (double)meta.col( test_pos[0] );
    dtest_pos[1] = (double)meta.row( test_pos[1] );
    cluster_info.closestHitsInCluster( clusterid, dtest_pos, src_data, meta, cm_per_tick, cm_per_wire, hitlist, 5 );
  }

  void StopMuTracker::getClosestHitsInList( const std::vector<int>& test_pos, const Hit2DList& src_list, const larcv::ImageMeta& meta, 
					    std::vector< std::pair<int,double> >& hitlist ) {
    // a wrapper around dbscanOutput::closestHitsInCluster
    // have to convert the 2D pixel col/row information into double test_point
    std::vector<double> dtest_pos(test_pos.size(),0);
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*m_config.usec_per_tick;
    float cm_per_wire = m_config.cm_per_wire;
    
    dtest_pos[0] = (double)meta.col( test_pos[0] );
    dtest_pos[1] = (double)meta.row( test_pos[1] );
    src_list.closestHits( dtest_pos, meta, cm_per_tick, cm_per_wire, hitlist, 5 );
  }

  // ------------------------------------------------------------------------------------------------
  // *** Re-orientation Routines ****
  // ------------------------------------------------------------------------------------------------

  // (write one)
  
  // ------------------------------------------------------------------------------------------------
  // *** UTILITY FUNCTIONS ****
  // ------------------------------------------------------------------------------------------------
  
  float StopMuTracker::_norm( std::vector<float>& vec ) {
    float norm = 0;
    for (size_t i=0; i<vec.size(); i++) norm += vec[i]*vec[i];
    norm = sqrt(norm);
    for (size_t i=0; i<vec.size(); i++) vec[i] /= norm;
    return norm;
  }

  void StopMuTracker::_wire2pixel( const int tick, const std::vector<int>& wid, const larcv::ImageMeta& meta, std::vector<int>& pixel_col, int& pixel_row ) {
    pixel_col.resize(m_config.nplanes,0.0);
    if ( pixel_col.size()!=wid.size() ) {
      throw std::runtime_error("StopMuTracker::_wire2pixel: provided wires for more or less than 3 planes");
    }
    for (int p=0; p<m_config.nplanes; p++)
      pixel_col[p] = meta.col( wid[p] );
    pixel_row = meta.row( tick );
  }

  std::vector<larcv::Image2D> StopMuTracker::fillSortedHit2Dlist( const larcv::ImageMeta& meta, 
								  const std::vector< std::vector<int> >& start2d,  const std::vector< std::vector<float> >& start_dir2d,
								  std::vector<Hit2DList>& hitlists, std::vector<int>& clusterid, const std::vector< float >& match_radius ) {
    // given a starting pixel, uses stored clusters (m_clusters) and list of dbPoints (m_imghits) 
    // to provide a list of Hit2D objects sorted by distance from start point
    // also returns an image containing which contains the hits that form the cluster in question

    std::cout << "StopMuTracker::fillSortedHit2Dlist" << std::endl;

    std::vector<larcv::Image2D> img_clusters;

    clusterid.resize(m_config.nplanes,-1);
    hitlists.clear();
    hitlists.resize(m_config.nplanes);
    for (int p=0; p<m_config.nplanes; p++) {
      std::vector<double> testpoint(2);
      testpoint[0] = start2d.at(p)[0]; // X
      testpoint[1] = start2d.at(p)[1]; // Y
      std::cout << "plane " << p << " testpoint: (col,row)=(" << testpoint[0] << "," << testpoint[1] << ")" 
		<< " (wire,tick)=(" << meta.pos_x( testpoint[0] ) << "," << meta.pos_y( testpoint[1] ) << ")" << std::endl;
      int match = m_clusters.at(p).findMatchingCluster( testpoint, m_imghits.at(p), match_radius[p] );
      std::cout << "plane " << p << " matching cluster=" << match << std::endl;
      clusterid[p] = match;
    }
    
    for (int p=0; p<m_config.nplanes; p++) {

      larcv::Image2D img_cluster( meta );
      img_cluster.paint(0.0);
      int match = clusterid.at(p);
      std::cout << "plane=" << p << " start_dir2d=(" << start_dir2d.at(p)[0] << "," << start_dir2d.at(p)[1] << ") "
		<< " size=" << m_clusters.at(p).clusters.at(match).size()
		<< std::endl;

      // we make a sorted list of pixels by distance
      // we also need an initial direction

      for (size_t ihit=0; ihit<m_clusters.at(p).clusters.at(match).size(); ihit++) {
	int hitidx = m_clusters.at(p).clusters.at(match).at(ihit);
	int x = m_imghits.at(p).at(hitidx)[0];
	int y = m_imghits.at(p).at(hitidx)[1];

	// dir from start to point. dir2d derives from cm-scales
	std::vector<float> dir(2);
	dir[0] =  (x-start2d.at(p)[0])*meta.pixel_width()*0.3; // cm
	dir[1] = -(y-start2d.at(p)[1])*meta.pixel_height()*m_config.usec_per_tick*::larutil::LArProperties::GetME()->DriftVelocity(); // cm
	float norm = sqrt( dir[0]*dir[0] + dir[1]*dir[1] );
	for (int i=0; i<2; i++) dir[i] /= norm;
	
	float cosine = 0.;
	for (int i=0; i<2; i++) cosine += dir[i]*start_dir2d.at(p)[i];

	Hit2D hit;
	hit[0] = x;
	hit[1] = y;
	if ( cosine>0 )
	  hit.distance = norm;
	else
	  hit.distance = -norm;
	if ( norm==0 ) {
	  hit.distance = 0;
	  cosine = 1;
	}

	if ( cosine>0 ) {
	  hitlists[p].emplace(std::move(hit));
	  hitlists[p].sort();
	}
	
	//std::cout << " hitidx=" << hitidx << ": (" << hit[0] << "," << hit[1] << ") cosine=" << cosine << std::endl;
	img_cluster.set_pixel((int)y,(int)x,250.0);
      }

      
      std::cout << "plane " << p << " number of hits=" << hitlists[p].size() << std::endl;
//       for (int ihit=0; ihit<(int)hitlists[p].size(); ihit++) {
// 	std::cout << " [#" << ihit << "] (" << hitlists[p].at(ihit)[0] << "," << hitlists[p].at(ihit)[1] << ")"
// 		  << " " << hitlists[p].at(ihit).distance << std::endl; 
//       }

      img_clusters.emplace_back( std::move(img_cluster) );
      
    }//end of loop over planes for sorting hits
    
    return img_clusters;
  }// StopMuTracker::fillSortedHit2Dlist

  // ------------------------------------------------------------------------------------------------
  // *** Hit2DList Functions ****
  // ------------------------------------------------------------------------------------------------

  void Hit2DList::closestHits( std::vector<double>& test_pos, const larcv::ImageMeta& meta, const float cm_per_tick, const float cm_per_wire,
			       std::vector< std::pair<int,double> >& hitlist, const int max_nhits, const int ignore_marked ) const {

    struct mycompare_t {
#if __GNUC__ >= 6
      bool operator() ( std::pair<int,double>& lhs, std::pair<int,double>& rhs ) { 
        if ( lhs.second<rhs.second ) 
          return true;
        return false;
      }
#else
      bool operator() ( std::pair<int,double> lhs, std::pair<int,double> rhs ) { 
        if ( lhs.second<rhs.second ) 
          return true;
        return false;
      }      
#endif
    } mycompare;

    hitlist.clear();
    for (int ihit=0; ihit<(int)size(); ihit++) {
      double hit_dist  = 0.0;
      const Hit2D& hitpos = at(ihit);
      if ( ignore_marked==1 && hitpos.marked==1 ) continue;
      double dt = ((double)hitpos[1]-test_pos[1])*meta.pixel_height()*cm_per_tick;
      double dw = ((double)hitpos[0]-test_pos[0])*meta.pixel_width()*cm_per_wire;
      hit_dist = sqrt( dt*dt + dw*dw ); // in cm
      std::pair<int,double> test_hit( ihit, hit_dist );
      if ( max_nhits>0 ) {
	// we have to care about the number of hits in the list
	if ( (int)hitlist.size()<max_nhits ) {
	  // but not now
	  hitlist.emplace_back( test_hit );
	}
	else {
	  bool isbetter = false;
	  for (size_t ilisthit=0; ilisthit<hitlist.size(); ilisthit++) {
	    if ( hitlist.at(ilisthit).second > hit_dist ) {
	      isbetter = true;
	      break;
	    }
	  }
	  if ( isbetter )
	    hitlist.emplace_back( test_hit );
	}
      }
      else {
	// dont care about the number of hits. add it to hit list.
	hitlist.emplace_back( test_hit );
      }
      
      // sort the current hitlist vector
      std::sort( hitlist.begin(), hitlist.end(), mycompare );
      // truncate the end
      if ( max_nhits>0 && max_nhits<(int)hitlist.size() )
	hitlist.resize(max_nhits);
            
    }//end of loop over hits

    if ( _verbose_ ) {
      std::cout << "sorted distances from test point (total hits=" << size() << "):" << std::endl;
      for (size_t i=0; i<hitlist.size(); i++) {
	  std::cout << " #" << i << ": " << hitlist.at(i).first << " " << hitlist.at(i).second << std::endl;
      }
    }

  }

  // ------------------------------------------------------------------------------------------------
  // *** String Minimizer ****
  // ------------------------------------------------------------------------------------------------

  FitDataHack* FitDataHack::_global_instance = NULL;

  double stop_mu_position_score(const double *xx )
  {
    // fit score comes from these components:
    // 1) how close it is to a pixel with charge above threshold
    // 2) how straight the step is with respect to the last step
    // 3) does the charge expected along the step consistent with that seen?
    
    // needs to find the closest hit
    Double_t cosz = xx[0];
    Double_t phi  = xx[1];
    bool _verbose_ = false;
    
    const larcv::ImageMeta& meta = *(FitDataHack::getMe()->data_meta);
    const std::vector<larcv::Image2D>& img_v = *(FitDataHack::getMe()->data_images);
    const std::vector<Hit2DList>& hitlists = *(FitDataHack::getMe()->data_hitlist);
    const std::vector<float>& anchor   = *(FitDataHack::getMe()->data_anchor);
    const std::vector<float>& prev_dir = *(FitDataHack::getMe()->data_prev_dir);
    const float dist_weight   = FitDataHack::getMe()->dist_weight;
    const float bend_weight   = FitDataHack::getMe()->bend_weight;
    const float charge_weight = FitDataHack::getMe()->charge_weight;
    const float step_size     = FitDataHack::getMe()->step_size_cm;

    // first make point and direction
    std::vector<float> dir(3,0.0);
    std::vector<float> pos(3,0.0);
    float rphi = sqrt(1-cosz*cosz);
    dir[2] = cosz;
    dir[1] = rphi*sin(phi);
    dir[0] = rphi*cos(phi);
    for (int i=0; i<3; i++) pos[i] = anchor[i] + step_size*dir[i];

    // get plane positions
    int tick;
    std::vector<int> wids;
    larlitecv::StopMuTracker::imagePositions( pos, tick, wids );
    if ( tick<meta.min_y() )  tick = meta.min_y();
    if ( tick>=meta.max_y() ) tick = meta.max_y()-1;
    
    int tick_anchor;
    std::vector<int> wids_anchor;
    larlitecv::StopMuTracker::imagePositions( anchor, tick_anchor, wids_anchor );
    if ( tick_anchor<meta.min_y() )  tick_anchor = meta.min_y();
    if ( tick_anchor>=meta.max_y() ) tick_anchor = meta.max_y()-1;

    // find the closest point on each plane
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    float cm_per_wire = 0.3; //< get this in here
    std::vector<float> closest_dists(3,1.0e3);
    float total_dist = 0.;
    if ( _verbose_ )
      std::cout << "test (" <<  cosz << "," << phi << ") distances from (" << pos[0] << "," << pos[1] << "," << pos[2] << "): ";
    for (int p=0; p<3; p++) {
      const Hit2DList& planehits = hitlists.at(p);
      std::vector< std::pair<int,double> > closest;
      std::vector<double> testpoint(2,0.0);
      testpoint[0] = meta.col( wids[p] );
      testpoint[1] = meta.row( tick );
      if ( testpoint[0]< 0 ) testpoint[0];
      if ( testpoint[0]>=meta.cols() ) testpoint[0] = meta.cols()-1;
      if ( testpoint[1]< 0 ) testpoint[1] = 0;
      if ( testpoint[1]>=meta.rows() ) testpoint[1] = meta.rows()-1;
      planehits.closestHits( testpoint, meta, cm_per_tick, cm_per_wire, closest, 3, 1 );
      if ( closest.size()>0 && closest.front().second<closest_dists[p]) {
	closest_dists[p] = closest.front().second;
	if ( _verbose_ )
	  std::cout << " p" << p << "=" << closest.front().second;
      }
      total_dist+= closest_dists[p]*closest_dists[p];
    }

    // cosine between previous and current direction
    float cosine = 0.;
    for (int i=0; i<3; i++) cosine += prev_dir[i]*dir[i];
    float cos_score = 1.0-cosine;
    if ( _verbose_ )
      std::cout << " bend cos=" << cosine;
    
    // charge along the path. 
    // we use Y-only because the calorimetry there is more reliable
    // what is the charge per y-pixel
    float mip_adc_per_voxel = 30.0;
    
    // first count the charge along the step
    int n_ypixels = 0;
    float tot_adc_on_yplane = 0.;
    std::vector<float> ydir(2);
    ydir[0] = float(meta.col(wids[2]))-float(meta.col(wids_anchor[2]));
    ydir[1] = float(meta.row(tick)) - float(meta.row(tick_anchor));
    float ynorm = sqrt( ydir[0]*ydir[0] + ydir[1]*ydir[1] );
    if ( fabs(ydir[0])<1.0e-3 && fabs(ydir[1])<1.0e-3 ) {
      int col = meta.col(wids[2]);
      int row = meta.row(tick);
      if ( col<0 ) col = 0;
      if ( col>=(int)meta.cols() ) col = (int)meta.cols()-1;
      if ( row<0 ) row = 0;
      if ( row>=(int)meta.rows() ) row = (int)meta.rows()-1;
      tot_adc_on_yplane = img_v[2].pixel( row, col );
      n_ypixels++;
    }
    else {
      // need to step through pixels
      //std::cout << " ynorm=" << ynorm << " ydir=(" << ydir[0] << "," << ydir[1] << ") dcol=" << meta.col(wids[2]) << "-" << meta.col(wids_anchor[2]) << " ";
      ydir[0] /= ynorm;
      ydir[1] /= ynorm;
      float ncols = fabs( float(meta.col(wids[2]))-float(meta.col(wids_anchor[2])) );
      float nrows = fabs( float(meta.row(tick)) - float(meta.row(tick_anchor)) );
      if ( _verbose_ )
	std::cout << " nrows=" << nrows << " ncols=" << ncols << std::endl;
      if ( fabs(ydir[0])>fabs(ydir[1]) ) {
	// cols faster moving than rows
	for (int icol=0; icol<=(int)ncols; icol++) {
	  int row = float(meta.row(tick_anchor)) + (ydir[1]/ydir[0])*icol;
	  int col = float(meta.col(wids_anchor[2])) + icol;
	  if ( col<0 ) col = 0;
	  if ( col>=(int)meta.cols() ) col = (int)meta.cols()-1;
	  if ( row<0 ) row = 0;
	  if ( row>=(int)meta.rows() ) row = (int)meta.rows()-1;
	  tot_adc_on_yplane += img_v[2].pixel( row, col );
	  n_ypixels++;
	}
      }
      else {
	// rows faster than cols
	for (int irow=0; irow<=(int)nrows; irow++) {
	  int col = float(meta.col(wids_anchor[2])) + (ydir[0]/ydir[1])*irow;
	  int row = float(meta.row(tick_anchor)) + irow;
	  if ( col<0 ) col = 0;
	  if ( col>=(int)meta.cols() ) col = (int)meta.cols()-1;
	  if ( row<0 ) row = 0;
	  if ( row>=(int)meta.rows() ) row = (int)meta.rows()-1;
	  tot_adc_on_yplane += img_v[2].pixel( row, col );
	  n_ypixels++;
	}	
      }
    }

    // now predict the charge for the step
    // basically, to do this, we need to create a voxel whose width is the wire pitch, the height is the dirft distance for one tick
    // and length is the wire. The latter is in principle, but in practice the track could end (since its stopping) any time 
    
    // ok, so we need to know which boundary the step hits first: (yz, zy, xz)
    float sx = (1.0*cm_per_tick)/dir[0];  // crossing over to the next tick
    float sy= (117.0-anchor[1])/dir[1]; // the actual cell
    float sz = (1.0*cm_per_wire)/dir[2];  // crossing over to the next wire
    if ( dir[1]<0 )
      sy = (-117.0-anchor[1])/dir[1]; 
    float s[3] = { sx, sy, sz };
    
    int shortest_dir = -1;
    for (int i=0; i<3; i++) {
      if ( dir[i]!=0 && (shortest_dir==-1 || shortest_dir>fabs(s[i]) ) ) {
	shortest_dir = i;
      }
    }
    
    // where does this point end up?
    std::vector<float> voxel_cross(3,0.0);
    float voxel_dist = 0.;
    for (int i=0; i<3; i++) {
      voxel_cross[i] = fabs(s[shortest_dir])*dir[i];
      voxel_dist += voxel_cross[i]*voxel_cross[i];
    }
    voxel_dist = sqrt(voxel_dist);
    
    float total_voxel_dist = voxel_dist*n_ypixels;
    
    float path_adc = mip_adc_per_voxel*(total_voxel_dist/cm_per_wire); // the ADC scale is assuming the wire crosses perpendicularly under the wire
    float adc_diff = path_adc-tot_adc_on_yplane;
    
    // final score
    if ( _verbose_ ) {
      std::cout << " voxel_dist=" << voxel_dist << " n_ypixels=" << n_ypixels;
      std::cout << " path_adc=" << path_adc << " sum_adc=" << tot_adc_on_yplane << " adc_diff=" << adc_diff << std::endl;
    }

    if ( _verbose_ )
      std::cout << std::endl;
    
    double score = 0.5*dist_weight*(total_dist) + 0.5*bend_weight*cos_score*cos_score + 0.5*charge_weight*adc_diff*adc_diff;
    
    return score;
  }

  void StopMuTracker::stopMuString( const std::vector<larcv::Image2D>& img_v, 
				    const std::vector< std::vector<int> >& start2d, const std::vector< std::vector<float> >& start_dir2d,
				    const std::vector< float >& start_pos3d, const std::vector<float>& start_dir3d, Step3D& trackstart ) {

    // parameters: move these to configuration file later
    float fStepSize_cm = 0.3; // 3 wires or ~10 microseconds or 20 ticks, which are a little more than 3 rows
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*m_config.usec_per_tick; // [cm/usec]*[usec/tick]
    float cm_per_wire = m_config.cm_per_wire;

    // find matching cluster
    std::vector<int> clusterid;
    const larcv::ImageMeta& meta = skel_v.at(0).meta();
    std::vector<Hit2DList> hitlists(3);
    std::vector<larcv::Image2D> img_clusters;
    std::vector<float> match_radius = {5,5,5};
    bool clustersfound = false;
    for (int itry=0; itry<3; itry++) {
      try {
	img_clusters = fillSortedHit2Dlist( meta, start2d, start_dir2d, hitlists, clusterid, match_radius );
      }
      catch (const std::exception& e) {
	std::cout << "could not find a cluster for the end point in each plane?" << std::endl;
      }
      clustersfound = true;
      for (int p=0; p<3; p++) {
	match_radius[p] *= 2.0;
	if ( clusterid[p]<0 ) {
	  clustersfound = false;
	}
      }
      if ( clustersfound ) break;
    }
    if ( !clustersfound ) {
      std::stringstream ss;
      ss << "no cluster error" << std::endl;
      throw std::runtime_error(ss.str());
    }

#ifdef USE_OPENCV
    //for debug output
    if ( m_verbosity>2 ) {
      for (int p=0; p<3; p++) {
    	const larcv::Image2D& img_cluster = img_clusters.at(p);
    	cv::Mat imgmat = larcv::as_mat( img_cluster );
    	std::stringstream ss;
    	ss << "baka_p" << p << ".jpg";
    	cv::imwrite( ss.str().c_str(), imgmat );
      }
    }
#endif

    // we need to check the quality of the cluster
    bool clusterok = true;
    for (size_t p=0; p<3; p++) {
      if ( hitlists.at(p).size()<1 ) {
	clusterok = false;
	break;
      }
    }

    if ( clusterok==false ) {
      return;
    }

    int istep = 0;
    bool isfinished = false;

    // we build 3D steps by fitting the position.  The step size is fixed, but the 3D angle is chosen. 
    // the potential is the closest distance to a hit + a "stress" term that wants to keep the line straight

    // current step, direction
    std::vector< float > current_pos =  start_pos3d;
    std::vector< float > current_dir =  start_dir3d;
    std::vector< float > proposed_pos(3,0.0);
    std::vector< float > proposed_dir(3,0.0);
    _norm( current_dir );

    // make the first Step
    trackstart.pos = current_pos;
    trackstart.dir = current_dir;
    trackstart.closesthits = start2d;
    trackstart.planepositions = start2d;

    // get the first step's plane position
    // (to do, not important at the moment)
    
    Step3D* current_step = &trackstart;

    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Scan");
    minimizer->SetMaxFunctionCalls(10000);
    minimizer->SetMaxIterations(10000); 
    minimizer->SetTolerance(0.001);
    minimizer->SetPrintLevel(0);
      
    ROOT::Math::Functor func(&stop_mu_position_score,2);
    double step[2] = {0.1,0.1};

    minimizer->SetFunction(func);

    while ( !isfinished && istep<1000 ) {
     
      if ( m_verbosity > 0 ) {
	std::cout << "===============================================================" << std::endl;
	std::cout << " Step " << istep << " [ptr=" << &(*current_step) << "]" << std::endl;
      }

      current_pos = current_step->pos;
      current_dir = current_step->dir;

      // update 3D step
      makeProposedPos( current_pos, current_dir, proposed_pos, fStepSize_cm );

      if ( m_verbosity > 0 ) {
	std::cout << " from current pos=(" << current_pos[0] << "," << current_pos[1] << "," << current_pos[2] << ")"
		  << " and current dir=(" << current_dir[0] << "," << current_dir[1] << "," << current_dir[2] << ")" << std::endl;
	std::cout << " stepsize=" << fStepSize_cm << " cm" << std::endl;
	std::cout << " made proposal: " << proposed_pos[0] << "," << proposed_pos[1] << "," << proposed_pos[2] << ")" << std::endl;
      }

      double min_value = 0.;
      double min_cosine = 2.0;
      if ( istep>=1 ) {
	// we need to minimize it
	
	// pass in the data we need
	FitDataHack::getMe()->data_images = &img_v;
	FitDataHack::getMe()->data_meta = &meta;
	FitDataHack::getMe()->data_hitlist = &hitlists;
	FitDataHack::getMe()->data_anchor = &current_pos;
	FitDataHack::getMe()->data_prev_dir = &( current_step->GetPrev().dir );
	FitDataHack::getMe()->dist_weight = m_config.fitter_distance_weight;
	FitDataHack::getMe()->bend_weight = m_config.fitter_bend_weight;
	FitDataHack::getMe()->charge_weight = 0.00;
	FitDataHack::getMe()->step_size_cm = m_config.fitter_step_size_cm;

	// set the initial values
	double variable[2] = { 0.98, 0.0 };
	// get the direction into cosz,phi
	variable[0] = current_dir[2];
	variable[1] = atan2(current_dir[1],current_dir[0]);

	minimizer->SetLimitedVariable(0, "cosz", variable[0], step[0], -1.0, 1.0 );
	minimizer->SetLimitedVariable(1,  "phi", variable[1], step[1], -3.14159, 3.14159 );

	// do the minimization
	if ( m_verbosity > 0 ) {
	  std::cout << "initialize minimizer with cosz=" << variable[0] << " phi=" << variable[1] << std::endl;
	  std::cout << "Run Minimizer." << std::endl;
	}
	minimizer->Minimize();

	const double *xs = minimizer->X();
	if ( m_verbosity > 0 )
	  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << minimizer->MinValue()  << std::endl;
	proposed_dir[2] = xs[0];
	proposed_dir[0] = sqrt(1.0-xs[0]*xs[0])*cos(xs[1]);
	proposed_dir[1] = sqrt(1.0-xs[0]*xs[0])*sin(xs[1]);
	_norm(proposed_dir);
	for (int i=0; i<3; i++) proposed_pos[i] = current_pos[i] + proposed_dir[i]*fStepSize_cm;
	min_value = minimizer->MinValue();
	if ( m_verbosity > 0 ) {
	  std::cout << "post-min pos: (" << proposed_pos[0] << "," << proposed_pos[1] << "," << proposed_pos[2] << ")" << std::endl;
	  std::cout << "post-min dir: (" << proposed_dir[0] << "," << proposed_dir[1] << "," << proposed_dir[2] << ")" << std::endl;
	  std::cout << "post-min func. value=" << min_value << std::endl;
	}
	min_cosine = 0.;
	for (int i=0; i<3; i++) min_cosine += current_dir[i]*proposed_dir[i];
      }
      else {
	// we skip the minimization for the first step in order to get a previous direction for input
	proposed_dir = current_dir;
      }

      // stopping condition test
      if ( m_verbosity > 0 )
	std::cout << "post-min cosine: " << min_cosine << std::endl;
      if ( min_cosine<-0.2 ) {
	if ( m_verbosity > 0 )
	  std::cout << "Stopping condition met." << std::endl;
	break;
      }

      // calculate 2D positions
      int tick;
      std::vector<int> wid;
      imagePositions( proposed_pos, tick, wid );
      if ( tick<meta.min_y() )  tick = meta.min_y();
      if ( tick>=meta.max_y() ) tick = meta.max_y()-1;

      if ( m_verbosity > 0 )
	std::cout << "proposal location on the 2D planes: tick=" << tick << " planes U=" << wid[0] << " V=" << wid[1] << " Y=" << wid[2] << std::endl;

      std::vector<int> pixel_cols;
      int pixel_row;
      _wire2pixel( tick, wid, meta, pixel_cols, pixel_row );
      if ( m_verbosity > 0 ) 
	std::cout << "proposal location on the image: row=" << pixel_row << " plane col U=" << pixel_cols[0] << " V=" << pixel_cols[1] << " Y=" << pixel_cols[2] << std::endl;


      // update tracker
      // prepare wire info
      if ( m_verbosity > 0 )
	std::cout << "update tracker/prepare wire info" << std::endl;
      std::vector<Point2D_t> closest_hits_list(3);
      std::vector<Point2D_t> plane_pos_list(3);
      bool close_enough = true;
      for (int p=0; p<3; p++) {
	std::vector<int> theplanehit(2);
	theplanehit[1] = pixel_row;
	theplanehit[0] = pixel_cols[p];
	plane_pos_list[p] = theplanehit;

	std::vector<double> dhit(2,0);
	dhit[0] = theplanehit[0];
	dhit[1] = theplanehit[1];
	std::vector<int> theclosesthit(2,0);
	std::vector< std::pair<int,double> > close_list;
	if ( m_verbosity>1 )
	  hitlists[p]._verbose_ = true;
	else
	  hitlists[p]._verbose_ = false;
	hitlists[p].closestHits( dhit, meta, cm_per_tick, cm_per_wire, close_list, 3, 1 );
	hitlists[p]._verbose_ = false;
	if ( m_verbosity > 0 )
	  std::cout << "mark hit up to: ";
	if ( close_list.size()>0 ) {
	  int hitidx = close_list.front().first;
	  hitlists[p].markUpTo( hitidx );
	  theclosesthit[0] = hitlists[p].at( hitidx )[0];
	  theclosesthit[1] = hitlists[p].at( hitidx )[1];
	  closest_hits_list[p] = theclosesthit;
	  if ( m_verbosity > 0 ) {
	    std::cout << " p" << p << "=" << hitidx 
		      << " closest hit=(" << theclosesthit[0] << "," << theclosesthit[1] << ")"
		      << std::endl;
	  }
	  if ( close_list.front().second>m_config.max_steppoint_dist_from_cluster ) {
	    close_enough = false;
	  }
	}
	else {
	  close_enough = false;
	}
      }
      
      Step3D* next_step = new Step3D( proposed_pos, proposed_dir, closest_hits_list, plane_pos_list, *current_step );
      current_step = next_step;
      if ( m_verbosity > 0 )
	std::cout << "updated Step3D list." << std::endl;

      bool hits_left = true;
      for (int p=0; p<3; p++) {
	if ( hitlists[p].back().marked==1 ) {
	  hits_left = false;
	  break;
	}
      }
	
      if ( !hits_left ) {
	if ( m_verbosity > 0 )
	  std::cout << "no more hits in the cluster." << std::endl;
	break;
      }
      if ( !close_enough ) {
	if ( m_verbosity > 0 )
	  std::cout << "we have lost our way." << std::endl;
	break;
      }

	
      istep++;
    }
    
    if ( m_verbosity>0 ) {
      std::cout << "[End of stopmu tracker. numebr of steps=" << istep << ".]" << std::endl;
      std::cout << "===============================================================" << std::endl;    
    }
  }
  
}
