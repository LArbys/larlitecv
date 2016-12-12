#include "StopMuAlgo.h"
#include <assert.h>

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv
#include "UBWireTool/UBWireTool.h"
#include "dbscan/DBSCANAlgo.h"


namespace larlitecv {

  StopMuAlgo::StopMuAlgo() : verbose(0) {}
  StopMuAlgo::~StopMuAlgo() {}

  void StopMuAlgo::runTimeStepTracker( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Pixel2D>& start, 
				       std::vector< std::vector<larcv::Pixel2D> >& pixellist, std::list< std::array<float,3> >& spacepoints ) {
    
    int rneighbor = 10;
    int cneighbor = 10;
    float fThreshold = 10;
    int neighborhood = 20;
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick]
    float ftrigger_time = 3200;

    // get starting directions near the start point
    //std::vector< std::vector<float> > start_dir2d;
    PlaneVec2D_t local_dir2d;
    Vec3D_t start_dir3d;
    Vec3D_t start_spacepoint;
    getStartDirection( img_v, badch_v, start, rneighbor, cneighbor, fThreshold, start_spacepoint, local_dir2d, start_dir3d );
    
    Vec3D_t local_dir3d;
    Vec2DList_t history_dir2d;
    Vec3DList_t history_dir3d;
    
    int idx_hist2d = 0; //< labels the pos to be filled
    int idx_hist3d = 0; //< labels the pos to be filled
    // we use these rotating indices instead of constantly removing and adding vector items
    
    int nsteps = 100;
    int nlocal_steps = 5;
    
    // set first step
    std::vector<larcv::Pixel2D> max_start_v;
    bool startok = findNeighborhoodMax( img_v, start, fThreshold, neighborhood, max_start_v );
    if ( !startok ) {
      std::cout << "could not start." << std::endl;
      return;
    }
    spacepoints.emplace_back( start_spacepoint );

    std::vector<larcv::Pixel2D> current = max_start_v; // a copy
    pixellist.push_back( current );

    if ( verbose>0 ) {
      std::cout << "=========[ START ]===================" << std::endl;
      std::cout << "2D Plane Pixels: U=(" << current.at(0).X() << "," << img_v.at(0).meta().pos_y(current.at(0).Y()) << ")"
		<< " V=(" << current.at(1).X() << "," << img_v.at(0).meta().pos_y(current.at(1).Y()) << ")"
		<< " Y=(" << current.at(2).X() << "," << img_v.at(0).meta().pos_y(current.at(2).Y()) << ")" << std::endl;
      std::cout << "local direction: " 
		<< "U=(" << local_dir2d.at(0)[0] << "," << local_dir2d.at(0)[1] << ") "
		<< "V=(" << local_dir2d.at(1)[0] << "," << local_dir2d.at(1)[1] << ") "
		<< "Y=(" << local_dir2d.at(2)[0] << "," << local_dir2d.at(2)[1] << ") " << std::endl;
    }

    for (int istep=0; istep<nsteps; istep++) {

      if (verbose>0)
	std::cout << "=========[STEP #" << istep << "]===================" << std::endl;
      // step method:
      // 1) from local_dir2d, calculate proposal for next step
      // 2) search near proposal for maximum pixel
      // 3) use it to define a hit
      // 4) update local_dir2d with new step
      // 5) update global_dir2d with new step
    
      // issues to address
      // multiple hits (due to intersecting tracks to a delta ray)
      // dead wires
      // horizontal portions of the track
      
      std::vector<larcv::Pixel2D> next_pix_v;
      bool found_proposal = getProposedStep( img_v, badch_v, current, 
					     fThreshold, neighborhood, next_pix_v, 
					     local_dir2d );
      
      if ( !found_proposal ) {
	std::cout << "proposal for next step could not be found." << std::endl;
	break;
      }
      
      // get spacepoint for proposal
      int crosses = 0;
      double triarea = 0;
      std::vector<float> poszy(2,00);
      std::vector<int> wids(3,0);
      for (int p=0;p<3;p++)
	wids[p] = img_v.at(p).meta().pos_x( next_pix_v.at(p).X() );
      larcv::UBWireTool::wireIntersection( wids, poszy, triarea, crosses );
      float posx = (img_v.at(0).meta().pos_y(next_pix_v.at(0).Y())-ftrigger_time)*cm_per_tick;

      Vec3D_t nextsp;
      nextsp[0] = posx;
      nextsp[1] = poszy[1];
      nextsp[2] = poszy[0];

      // make 3D step direction
      Vec3D_t lastdir3d;
      const Vec3D_t& last_sp = spacepoints.back();
      float lastdir3d_norm = 0.;
      for (int i=0; i<3; i++) {
	lastdir3d[i] = nextsp[i] - last_sp[i];
	lastdir3d_norm += lastdir3d[i]*lastdir3d[i];
      }
      lastdir3d_norm = sqrt(lastdir3d_norm);
      for (int i=0; i<3; i++)
	lastdir3d[i] /= lastdir3d_norm;

      // check if the new step is going backwards
      if ( history_dir3d.size()>1 && amGoingBackwards( nextsp, history_dir3d, spacepoints, lastdir3d, -0.8, local_dir3d )) {
	if ( verbose> 0 ) std::cout << "going backwards." << std::endl;
	break;
      }

      // everything OK: make updates to position and directions
      
      // update 2D direction
      PlaneVec2D_t last_dir2d;
      for (int p=0; p<3; p++) {
	Vec2D_t pdir2d;
	pdir2d[0] = (float)next_pix_v.at(p).X() - (float)current.at(p).X();
	pdir2d[1] = (float)next_pix_v.at(p).Y() - (float)current.at(p).Y();
	float pdir2d_norm = sqrt( pdir2d[0]*pdir2d[0] + pdir2d[1]*pdir2d[1] );
	for (int i=0; i<2; i++ )
	  pdir2d[i] /= pdir2d_norm;
	//std::cout << "p=" << p << " adding " << pdir2d[0] << "," << pdir2d[1] << std::endl;
	last_dir2d.emplace_back( pdir2d );
      }
      
      // update the dir2d history
      history_dir2d.emplace_back( last_dir2d );
      if ( history_dir2d.size()>nlocal_steps )
	history_dir2d.pop_front(); // remove the first and oldest element

      if ( history_dir2d.size()>=nlocal_steps ) {
	for ( int p=0; p<3; p++ ) {
	  for (int i=0; i<2; i++)
	    local_dir2d.at(p)[i] = 0.0;
	}
	// use history to update the local dir2d
	for (int p=0; p<3; p++) {
	  for ( auto const& old_dir2d : history_dir2d ) {
	    for (int i=0; i<2; i++) {
	      local_dir2d.at(p)[i] += old_dir2d.at(p)[i]/float(history_dir2d.size());
	    }
	    //std::cout << "summing p=" << p << ": " << old_dir2d.at(p)[0] << "," << old_dir2d.at(p)[1] << std::endl;
	  }
	  //std::cout << "local_dir2d sum=" << local_dir2d.at(p)[0] << "," << local_dir2d.at(p)[1] << std::endl;
	}
	// normalize
	for (int p=0; p<3; p++) {
	  float local_dir2d_norm = 0.;
	  for (int i=0; i<2; i++) 
	    local_dir2d_norm += local_dir2d.at(p)[i]*local_dir2d.at(p)[i];
	  local_dir2d_norm = sqrt(local_dir2d_norm);
	  for (int i=0; i<2; i++)
	    local_dir2d.at(p)[i] /= local_dir2d_norm;
	}// whew
      }

      // update 3D pos
      spacepoints.emplace_back( nextsp );
      // update 3D direction
      history_dir3d.emplace_back( lastdir3d );
      if ( history_dir3d.size()>nlocal_steps )
	history_dir3d.pop_front();
      aveVec3DList( history_dir3d, local_dir3d );
      
      // set current to next point
      current.swap( next_pix_v );
      pixellist.push_back( current );

      // output

      // 2D
      if ( verbose>0 ) {
	std::cout << "2D Plane Pixels: U=(" << current.at(0).X() << "," << img_v.at(0).meta().pos_y(current.at(0).Y()) << ")"
		  << " V=(" << current.at(1).X() << "," << img_v.at(0).meta().pos_y(current.at(1).Y()) << ")"
		  << " Y=(" << current.at(2).X() << "," << img_v.at(0).meta().pos_y(current.at(2).Y()) << ")" << std::endl;
	std::cout << "local direction: " 
		  << "U=(" << local_dir2d.at(0)[0] << "," << local_dir2d.at(0)[1] << ") "
		  << "V=(" << local_dir2d.at(1)[0] << "," << local_dir2d.at(1)[1] << ") "
		  << "Y=(" << local_dir2d.at(2)[0] << "," << local_dir2d.at(2)[1] << ") " << std::endl;
	
	// 3D
	const Vec3D_t& newsp = spacepoints.back();
	const Vec3D_t& newdir3d = history_dir3d.back();
	std::cout << "3D pos=(" << newsp[0] << "," << newsp[1] << "," << newsp[2] << ")  triarea=" << triarea << std::endl;
	std::cout << "3D last dir=(" << newdir3d[0] << "," << newdir3d[1] << "," << newdir3d[2] << ")" << std::endl;
	std::cout << "3D local dir=(" << local_dir3d[0] << "," << local_dir3d[1] << "," << local_dir3d[2] << ")" << std::endl;
      }

      // go to next step
    }
    
  }
  
  bool StopMuAlgo::getProposedStep( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Pixel2D>& current,
				    float threshold, int neighborhood, std::vector<larcv::Pixel2D>& next_pix_v,
				    std::vector< std::array<float,2> >& local_dir2d ) {
    std::vector< larcv::Pixel2D > proposed_origin( current );
    bool found_proposal = false;
    for (int itry=0; itry<100; itry++) {
      if ( itry==50 ) {
	// going to flip time direction to see if it helps
	for (size_t p=0; p<local_dir2d.size(); p++) {
	  local_dir2d[p][1] *= -1.0;
	}
      }
      std::vector<larcv::Pixel2D> proposed_pix_v;
      bool ok = calcNextPos( img_v, local_dir2d, proposed_origin, proposed_pix_v );
      std::vector<larcv::Pixel2D> max_pix_v;
      ok = findNeighborhoodMax( img_v, proposed_pix_v, threshold, neighborhood, max_pix_v );
      if (!ok) {
	for (int ihit=0; ihit<max_pix_v.size(); ihit++) {
	  const larcv::Pixel2D& maxhit = max_pix_v.at(ihit);
	  if ( maxhit.Intensity()>=0 ) continue; // this pixel is OK
	  int dcol = 1.0;
	  if ( itry >= 50 ) dcol *= -1; // going backwards
	  if ( local_dir2d.at(ihit)[0]<0 ) dcol *= 1; // go in reverse
	  int adjusted_x = proposed_origin.at(ihit).X()+dcol;
	  int adjusted_y = proposed_origin.at(ihit).Y();
	  larcv::Pixel2D adjusted_pix( adjusted_x, adjusted_y );
	  adjusted_pix.Intensity(proposed_origin.at(ihit).Intensity());
	  // replace
	  proposed_origin[ihit] = adjusted_pix;
	}//end of loop over max hits
	continue; // go to next try
      }
      // if neighborhood maximum found
      std::vector< HitNeighborhood > hits_v;
      defineNeighborHoodHitSimple( img_v, max_pix_v, threshold, hits_v );
      
      // convert to pixels
      next_pix_v.clear();
      for ( auto& hits : hits_v ) {
	larcv::Pixel2D pix( hits.maxcol, max_pix_v.at(0).Y() );
	next_pix_v.emplace_back(pix);
      }
      // found a set of proposal pixels
      found_proposal = true;
      break;
    }
    return found_proposal;
  }

  void StopMuAlgo::getStartDirection( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
				      const std::vector<larcv::Pixel2D>& start, int rneighbor, int cneighbor, float fThreshold,
				      std::array<float,3>& start_spacepoint, std::vector< std::array<float,2> >& start_dir2d, std::array<float,3>& start_dir3d  ) {
				      
    
    float ftrigger_tick  = 3200.0;
    float step_size = 3.0;
    float max_time_step_cm = 1.0; // cm    
    
    const int nplanes = img_v.size();
    std::vector<dbscan::dbPoints> combo_points(nplanes);
    std::vector<dbscan::dbscanOutput> clustering_out;
    int endpoint_hitindex[nplanes];
    
    for (int p=0; p<nplanes; p++) {
      const larcv::ImageMeta& meta = img_v.at(p).meta();
      const larcv::Pixel2D& endpoint = start.at(p);
      //std::cout << "plane=" << p << " start point: (" << endpoint.X() << "," << endpoint.Y() << ")" << std::endl;
      // gather pixels above threshold that are around endpoint
      for (int dr=-rneighbor;dr<=rneighbor; dr++) {
	int row = endpoint.Y() + dr;
	if ( row<0 || row >= meta.rows() ) continue;
	for ( int dc=-cneighbor;dc<=cneighbor; dc++ ) {
	  int col = endpoint.X() + dc;
	  if ( col<0 || col>= meta.cols() ) continue;
	  if ( img_v.at(p).pixel( row, col ) > fThreshold || (dc==0 && dr==0) ) {
	    // if pixel above threshold
	    std::vector<double> point(2);
	    point[0] = (double)col;
	    point[1] = (double)row;
	    combo_points[p].emplace_back( point );
	    if ( dc==0 && dr==0 ) {
	      endpoint_hitindex[p] = combo_points[p].size()-1;
	    }
	  }
	}
      }
      // cluster points
      dbscan::DBSCANAlgo algo;
      dbscan::dbscanOutput clout = algo.scan( combo_points[p], 5, 3 );
      //std::cout << "neighborhood on plane=" << p <<  " has " <<  clout.clusters.size() << " clusters" << std::endl;
      // store clustering results
      clustering_out.emplace_back( clout );
    }//end of plane loop
    

    //  (b) we get the cluster that the point is attached to
    int pointclusters[nplanes];
    for (int p=0; p<nplanes; p++) {
      pointclusters[p] = -1;
      for (int c=0; c<clustering_out.at(p).clusters.size(); c++) {
	if ( pointclusters[p]!=-1 ) break; // already found the cluster, move on
	const std::vector<int>& cluster = clustering_out.at(p).clusters.at(c);
	for (int ihit=0; ihit<cluster.size(); ihit++) {
	  int hitidx = cluster.at(ihit);
	  if ( hitidx==endpoint_hitindex[p] ) {
	    pointclusters[p] = c;
	    break;
	  }
	}
      }
      //std::cout << "end point for plane=" << p << " is on cluster=" << pointclusters[p]
      //	<< " with " << clustering_out.at(p).clusters.at(pointclusters[p]).size() << " hits" << std::endl;
    }//end of loop over planes


    //  (c) we find the furthest pixel and get the direction
    float plane_dir[nplanes][2];
    float maxdist[nplanes];
    for (int p=0; p<nplanes; p++) {
      maxdist[p] = 0.0;
      const larcv::ImageMeta& meta = img_v.at(p).meta();
      plane_dir[p][0] = plane_dir[p][1] = 0;
      const std::vector<int>& cluster = clustering_out.at(p).clusters.at(pointclusters[p]);
      float hitpos[2] = {0};
      for (int ihit=0; ihit<cluster.size(); ihit++) {
	float dx = (combo_points[p].at(ihit)[0]-start.at(p).X())*meta.pixel_width()*0.3; // cm
	float dy = -(combo_points[p].at(ihit)[1]-start.at(p).Y())*meta.pixel_height()*0.5*::larutil::LArProperties::GetME()->DriftVelocity(); // cm 
	// (note, above negative, because of reverse time order used in image 2D)
	float dist = sqrt( dx*dx + dy*dy );
	if ( dist>maxdist[p] ) {
	  maxdist[p] = dist;
	  hitpos[0] = combo_points[p].at(ihit)[0];
	  hitpos[1] = combo_points[p].at(ihit)[1];
	  plane_dir[p][0] = dx/dist;
	  plane_dir[p][1] = dy/dist;
	}
      }
      //std::cout << "plane=" << p << " cluster dir=(" << plane_dir[p][0] << "," << plane_dir[p][1] << ") maxdist=" << maxdist[p] << " cm" << std::endl;
      std::array<float,2> sdir;
      sdir[0] = plane_dir[p][0];
      sdir[1] = plane_dir[p][1];
      start_dir2d.emplace_back( sdir );
    }//end of loop over planes


    //  (d) we use the 2D directions to provide a space point
    // find initial distance to move one tick, use that as guide to find target tick, picking the shortest one.
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick]
    float init_dist[nplanes];
    float target_tick[nplanes];
    float closest_tick = 1e6;
    for (int p=0; p<nplanes; p++) {
      if ( plane_dir[p][1]!=0.0 )
	init_dist[p] = fabs((step_size*img_v.at(p).meta().pixel_height()*cm_per_tick)/plane_dir[p][1]);// how many unit vectors required to go the step size
      else
	init_dist[p] = max_time_step_cm/cm_per_tick; 
      //std::cout << "init dist[" << p << "] " << init_dist[p] << " vs. " << max_time_step_cm/cm_per_tick << std::endl;
      if ( init_dist[p]>max_time_step_cm/cm_per_tick )
	init_dist[p] = max_time_step_cm/cm_per_tick;
      target_tick[p] = img_v.at(p).meta().pos_y( start.at(p).Y() ) + init_dist[p]*(plane_dir[p][1]/cm_per_tick);
      if ( target_tick[p]<closest_tick )
	closest_tick = target_tick[p];
      //std::cout << "initial direction determined by going from tick=" << img_v.at(p).meta().pos_y(start.at(p).Y()) << " --> " << closest_tick << std::endl;
    }// end of plane loop

    // get the wire coordinate of the target tick
    std::vector<int> target_col(nplanes,-1);
    std::vector<int> target_wires(nplanes,-1);
    std::vector<int> start_wires(nplanes,-1);
    for (int p=0; p<nplanes; p++) {
      if ( plane_dir[p][1]!=0.0 ) {
	float start_tick = img_v.at(p).meta().pos_y( start.at(p).Y() );
	//float s = (target_tick[p]-start_tick)*cm_per_tick/plane_dir[p][1];
	float s = (closest_tick-start_tick)*cm_per_tick/plane_dir[p][1];
	float dcm = plane_dir[p][0]*s; // change in cm's in wire direction after 's' units
	float dwire  = dcm/0.3; // [cm] / [cm/wire]
	float dcol   = dwire/img_v.at(p).meta().pixel_width();  // [dwire] / [dwire/pixel col]
	target_col[p] = (int)(start.at(p).X() + dcol);
      }
      else {
	target_col[p] = start.at(p).X();
      }
      start_wires[p]  = img_v.at(p).meta().pos_x( start.at(p).X() );
      target_wires[p] = img_v.at(p).meta().pos_x( target_col[p] );
    }
    //     std::cout << "initial column directions: " 
    // 	      << " u: " << start.at(0).X() << " --> " << target_col[0]
    // 	      << " v: " << start.at(1).X() << " --> " << target_col[1]
    // 	      << " y: " << start.at(2).X() << " --> " << target_col[2]
    // 	      << std::endl;
    
    // get the 3D position
    std::vector<float> spacepoint[2];
    double triarea[2] = {0};
    int crosses[2] = {0};
    // note these return (z,y,x), because, of course they do.
    larcv::UBWireTool::wireIntersection( start_wires,  spacepoint[0], triarea[0], crosses[0] );
    larcv::UBWireTool::wireIntersection( target_wires, spacepoint[1], triarea[1], crosses[1] );

    // use 3d points to get initial direction and position
    float start_tick = img_v.at(0).meta().pos_y( start.at(0).Y() );
    std::vector<float> init_pos(3,0.0);
    std::vector<float> end_pos(3,0.0);
    std::vector<float> init_dir(3,0.0);
    init_pos[0] = (start_tick-ftrigger_tick)*cm_per_tick;
    init_pos[1] = spacepoint[0][1];
    init_pos[2] = spacepoint[0][0];
    end_pos[0]  = (closest_tick-ftrigger_tick)*cm_per_tick;
    end_pos[1]  = spacepoint[1][1];
    end_pos[2]  = spacepoint[1][0];
    //std::cout << "start 3D point: (" << init_pos[0] << "," << init_pos[1] << "," << init_pos[2] << ") area=" << triarea[0] << "(" << crosses[0] << ")" << std::endl;
    //std::cout << "end 3D point:   (" << end_pos[0] << "," << end_pos[1] << "," << end_pos[2] << ") area=" << triarea[1] << "(" << crosses[1] << ")" << std::endl;

    for (int i=0; i<3; i++) {
      start_spacepoint[i] = init_pos[i];
    }
    
    float start_norm = 0.0;
    for (int i=0; i<3; i++) {
      init_dir[i] = end_pos[i]-init_pos[i];
      start_norm += init_dir[i]*init_dir[i];
    }
    start_norm = sqrt(start_norm);
    for (int i=0; i<3; i++) { 
      init_dir[i] /= start_norm;
      start_dir3d[i] = init_dir[i];
    }

  }//end of getStartDirection

  void StopMuAlgo::getStartDirectionV( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
				       const std::vector<larcv::Pixel2D>& start, int rneighbor, int cneighbor, float fThreshold,
				       std::vector<float>& start_spacepoint, std::vector< std::vector<float> >& start_dir2d, std::vector<float>& start_dir3d ) {
    // wraps the above function, because python bindings can't handle array yet
    Vec3D_t start_sp;
    PlaneVec2D_t start2d;
    Vec3D_t start3d;
    getStartDirection( img_v, badch_v, start, rneighbor, cneighbor, fThreshold, start_sp, start2d, start3d );
    start_spacepoint.resize(3,0.0);
    start_dir3d.resize(3,0.0);
    for (int i=0; i<3; i++) {
      start_spacepoint[i] = start_sp[i];
      start_dir3d[i] = start3d[i];
    }

    for (int p=0; p<3; p++) {
      std::vector<float> dir2d(2,0.0);
      for (int i=0; i<2; i++) {
	dir2d[i] = start2d.at(p)[i];
      }
      start_dir2d.emplace_back( dir2d );
    }
  }

  void StopMuAlgo::defineNeighborHoodHitSimple( const std::vector<larcv::Image2D >& img_v, const std::vector<larcv::Pixel2D>& start_pix_v, float threshold,
						std::vector< HitNeighborhood >& hits ) {
    
    for (int p=0; p<img_v.size(); p++) {
      const larcv::ImageMeta& meta = img_v.at(p).meta();
      int dc = 0;
      int row = start_pix_v.at(p).Y();
      float maxval = 0;
      int maxcol = 0;
      while (true) {
	int col = start_pix_v.at(p).X()+dc;
	if (col<0 || col>=meta.cols()) break;
	float q = img_v.at(p).pixel(row,col);
	if (q>maxval) {
	  maxcol = col;
	  maxval = q;
	}
	if (q<threshold) {
	  // maybe should fall below a few times
	  break;
	}
	dc+=-1;
      }
      int start = start_pix_v.at(p).X()+dc;
      dc = 0;
      while (true) {
	int col = start_pix_v.at(p).X()+dc;
	if (col<0 || col>=meta.cols()) break;
	float q = img_v.at(p).pixel(row,col);
	if (q>maxval) {
	  maxcol = col;
	  maxval = q;
	}
	if (q<threshold) {
	  break;
	}
	dc+=+1;
      }
      int end = start_pix_v.at(p).X()+dc;
      HitNeighborhood hit( start, end, maxcol, maxval );
      hits.emplace_back(hit);
    }//end of loop over planes
  }//end of loop over 'defineNeighborHoodHitSimple'
  
  bool StopMuAlgo::calcNextPos( const std::vector<larcv::Image2D>& img_v, const std::vector< std::array<float,2> >& dir2d, const std::vector< larcv::Pixel2D >& pix_v,
				std::vector<larcv::Pixel2D>& next_pix_v ) {
    for (int p=0; p<img_v.size(); p++) {
      if ( dir2d.at(p)[1]==0 ) {
	std::cout << "no time component for dir" << std::endl;
	// this is a horizontal track
	assert(false);
	return false;
      }
      int  rowsign = 1;
      if (dir2d.at(p)[1]<0)
	rowsign = -1;
      const larcv::Pixel2D& hit = pix_v.at(p);
      int row = hit.Y();
      int hit_row = row + rowsign;
      int col = int(hit.X() + dir2d.at(p)[0]*(1.0/dir2d.at(p)[1]));
      if ( col<0 ) 
	col = 0;
      if ( col>=img_v.at(p).meta().cols() ) 
	col = img_v.at(p).meta().cols()-1;
      larcv::Pixel2D pix( col, hit_row );
      next_pix_v.emplace_back(pix);
    }
    return true;
  }

  bool StopMuAlgo::findNeighborhoodMax( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Pixel2D>& pix_v, float threshold, int neighborhood,
					std::vector<larcv::Pixel2D>& max_pix_v ) {
    bool ok = true;
    for (int p=0; p<img_v.size(); p++) {
      int col = pix_v.at(p).X();
      int row = pix_v.at(p).Y();
      float maxval = -1.0;
      float maxcol = 0;
      for (int c=col-neighborhood; c<=col+neighborhood; c++) {
	if ( c<0 || c>=img_v.at(p).meta().cols() ) 
	  continue;
	float q = img_v.at(p).pixel(row,c);
	if (q<threshold) 
	  continue;
	if (q>maxval) {
	  maxval = q;
	  maxcol = c;
	}
      }
      //print "plane=",p," maxcol at row=",row," col=",c," maxval=",maxval
      larcv::Pixel2D maxpix( (unsigned int)maxcol, row );
      maxpix.Intensity(maxval);
      max_pix_v.emplace_back( maxpix );
      if ( maxval<0 )
	ok = false;
    }
    return ok;
  }

  bool StopMuAlgo::findNearestPix( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Pixel2D>& pix_v, const PlaneVec2D_t& dir2d,
				   const float threshold, const int neighborhood,
				   std::vector<larcv::Pixel2D>& max_pix_v ) {
    /// fix me
    bool ok = true;
    for (int p=0; p<img_v.size(); p++) {
      int col = pix_v.at(p).X();
      int row = pix_v.at(p).Y();
      float maxval = -1.0;
      float maxcol = 0;
      for (int c=col-neighborhood; c<=col+neighborhood; c++) {
	if ( c<0 || c>=img_v.at(p).meta().cols() ) 
	  continue;
	float q = img_v.at(p).pixel(row,c);
	if (q<threshold) 
	  continue;
	if (q>maxval) {
	  maxval = q;
	  maxcol = c;
	}
      }
      //print "plane=",p," maxcol at row=",row," col=",c," maxval=",maxval
      larcv::Pixel2D maxpix( (unsigned int)maxcol, row );
      maxpix.Intensity(maxval);
      max_pix_v.emplace_back( maxpix );
      if ( maxval<0 )
	ok = false;
    }
    return ok;
  }

  bool StopMuAlgo::amGoingBackwards( const Vec3D_t& nextsp, const Vec3DList_t& history_dir3d, const Vec3DList_t& spacepoints, const Vec3D_t& lastdir3d, 
				     const float fbackward_cos_threshold, Vec3D_t& local_dir3d ) {
				 
    // determines if going backward 
    // also calculates lastdir3d and local_dir3d
    const Vec3D_t& lastsp = spacepoints.back();
    
    // calculate average of vectors in history_dir3d
    aveVec3DList( history_dir3d, local_dir3d );
    
    float cos3d = 0.0;
    for (int i=0; i<3; i++)
      cos3d += lastdir3d[i]*local_dir3d[i];
    if ( verbose>0 )
      std::cout <<  "checking 3D cosine: " << cos3d << std::endl;
    if (cos3d<fbackward_cos_threshold)
      return true;
    return false;
  }

  void StopMuAlgo::aveVec3DList( const Vec3DList_t& history_dir3d, Vec3D_t& ave_dir3d ) {
    // calculate average of vectors in history_dir3d
    for (int i=0; i<3; i++)
      ave_dir3d[i] = 0.0;
    for (auto const& old_dir3d : history_dir3d ) {
      for (int i=0; i<3; i++)
	ave_dir3d[i] += old_dir3d[i]/float(history_dir3d.size());
      float norm2 = 0.;
      for (int i=0; i<3; i++)
	norm2 += ave_dir3d[i]*ave_dir3d[i];
      norm2 = sqrt(norm2);
      for (int i=0; i<3; i++)
	ave_dir3d[i] /= norm2;
    }
    
  }

}
