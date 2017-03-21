#include "AStarNodes2BMTrackCluster3D.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

namespace larlitecv {
	
  BMTrackCluster3D AStarNodes2BMTrackCluster3D( const std::vector<AStar3DNode>& path, 
    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
    const std::vector<const larcv::Pixel2D*>& start_pt, const std::vector<const larcv::Pixel2D*>& end_pt,
    const int pixel_tag_neighborhood, const float link_step_size, const std::vector<float>& pixel_thresholds ) {

    const larcv::ImageMeta& meta = img_v.front().meta();
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick] 
    const int nplanes = img_v.size();

    // fill track3d data
    BMTrackCluster3D track3d;

    // Start Point Information
    track3d.start_type = (larlitecv::BoundaryEnd_t) int(start_pt.front()->Intensity());
    track3d.row_start  = start_pt.front()->Y();
    track3d.tick_start = img_v.front().meta().pos_y( track3d.row_start );
    track3d.start_wire.resize(nplanes,0);
    track3d.start3D.resize(nplanes,0);
    for (int i=0; i<nplanes; i++) {
      track3d.start3D[i] = path.back().tyz[i];
      track3d.start_wire[i] = meta.pos_x( start_pt[i]->X() );
    }   
    track3d.start3D[0] = (track3d.start3D[0]-3200)*cm_per_tick;

    // End Point Information
    if ( end_pt.size()!=img_v.size() ) {
    	// No end specified
      track3d.end_type   = larlitecv::kUndefined;
      track3d.tick_end   = path.front().tyz.at(0);
      track3d.row_end    = meta.row( track3d.tick_end );
      track3d.end_wire.resize(nplanes,0);
      track3d.end3D.resize(nplanes,0);
      track3d.end3D[0] = (track3d.end3D[0]-3200)*cm_per_tick;    
      for (int i=1; i<nplanes; i++)
        track3d.end3D[i] = path.front().tyz[i];
      Double_t xyz_end[3];
      for (int i=0; i<nplanes; i++) 
        xyz_end[i] = track3d.end3D[i];    
      for (int i=0; i<nplanes; i++) {
        float fwire = larutil::Geometry::GetME()->WireCoordinate(xyz_end,i);
        fwire = ( fwire<0 ) ? 0 : fwire;
        fwire = ( fwire>=meta.max_x() ) ? meta.max_x()-1.0 : fwire;
        track3d.end_wire[i] = (int)fwire;
      }
    }
    else {
      track3d.start_type = (larlitecv::BoundaryEnd_t) int(end_pt.front()->Intensity());
      track3d.row_end  = end_pt.front()->Y();
      track3d.tick_end = img_v.front().meta().pos_y( track3d.row_end );
      track3d.end_wire.resize(nplanes,0);
      track3d.end3D.resize(nplanes,0);
      for (int i=0; i<nplanes; i++) {
        track3d.end3D[i] = path.back().tyz[i];
        track3d.end_wire[i] = meta.pos_x( end_pt[i]->X() );
      }   
      track3d.end3D[0] = (track3d.end3D[0]-3200)*cm_per_tick;    	
    }


    // Prepare Track2D objects and an empty image to track which pixels we've marked
    std::vector<larcv::Image2D> tagged_v;
    for (int p=0; p<nplanes; p++) {
      BMTrackCluster2D track2d;
      BoundaryEndPt start( track3d.row_start, meta.col( track3d.start_wire[p]), track3d.start_type );
      track2d.start = start;

      if ( end_pt.size()!=img_v.size() ) {
        BoundaryEndPt end( track3d.row_end, meta.col( track3d.end_wire[p] ), larlitecv::kUndefined );
        track2d.end   = end;        
      }
      else {
        BoundaryEndPt end( track3d.row_end, meta.col( track3d.end_wire[p]), track3d.end_type );
        track2d.end = end;
      }

      track2d.plane = p;
      track3d.plane_paths.emplace_back( std::move(track2d) );
      larcv::Image2D tagged( img_v.at(p).meta() );
      tagged.paint(0.0);
      tagged_v.emplace_back( std::move(tagged) );
    }

    float nbad_nodes = 0;
    float total_nodes = 0;
    int nnodes = (int)path.size();
    for ( int inode=nnodes-1; inode>=1; inode-- ) {

      const AStar3DNode& node      = path.at(inode);
      const AStar3DNode& next_node = path.at(inode-1);
      if ( node.badchnode )
        nbad_nodes+=1.0;
      total_nodes+=1.0;

      float dir3d[3];
      float step0[3];
      float dist = 0.;
      for (int i=0; i<3; i++) {
        dir3d[i] = next_node.tyz[i] - node.tyz[i];
        step0[i] = node.tyz[i];
      }
      dir3d[0] *= cm_per_tick;
      step0[0] = (step0[0]-3200.0)*cm_per_tick;
      for (int i=0; i<3; i++)
        dist += dir3d[i]*dir3d[i];
      dist = sqrt(dist);
      for (int i=0; i<3; i++)
        dir3d[i] /= dist;

      int nsteps = dist/link_step_size+1;
      float stepsize = dist/float(nsteps);

      for (int istep=0; istep<=nsteps; istep++) {
        Double_t xyz[3];
        std::vector<double> pt(3,0.0);
        for (int i=0; i<3; i++) {
          xyz[i] = step0[i] + stepsize*istep*dir3d[i];
          pt[i] = xyz[i];
        }
        float tick = xyz[0]/cm_per_tick + 3200.0;
        if ( tick<=meta.min_y() || tick>=meta.max_y() ) continue;
        track3d.path3d.emplace_back( std::move(pt) );        
        int row = meta.row( tick );        
        std::vector<int> cols(3);
        for (int p=0; p<nplanes; p++) {
          float fwire = larutil::Geometry::GetME()->WireCoordinate(xyz,p);
          fwire = ( fwire<0 ) ? 0 : fwire;
          fwire = ( fwire>=meta.max_x() ) ? meta.max_x()-1.0 : fwire;
          cols[p] = meta.col( fwire );

          for (int dr=-pixel_tag_neighborhood; dr<=pixel_tag_neighborhood; dr++) {
            int r = row+dr;
            if ( r<0 || r>=(int)meta.rows()) continue;
            for (int dc=-pixel_tag_neighborhood; dc<=pixel_tag_neighborhood; dc++) {
              int c = cols[p]+dc;
              if ( c<0 || c>=(int)meta.cols()) continue;
              // tag pixels that are (1) untagged && (2) above threshold or bad channels
              if ( tagged_v.at(p).pixel(r,c)==0 && (img_v.at(p).pixel(r,c)>pixel_thresholds.at(p) || badch_v.at(p).pixel(r,c)>0 ) ) {
                tagged_v.at(p).set_pixel(r,c,255);
                larcv::Pixel2D pix(c,r);
                track3d.plane_paths.at(p).pixelpath.emplace_back( std::move(pix) );
              }
            }
          }//end of row loop
        }//end of plane loop
      }//end of steps

    }//end of node loop

    return track3d;
  }

}
