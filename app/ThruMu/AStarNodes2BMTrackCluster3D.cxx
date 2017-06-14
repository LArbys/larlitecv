#include "AStarNodes2BMTrackCluster3D.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv (crb)
#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  BMTrackCluster3D AStarNodes2BMTrackCluster3D( const std::vector<AStar3DNode>& path, const std::vector<larcv::Image2D>& img_v,
    const BoundarySpacePoint& start_pt, const BoundarySpacePoint& end_pt, const float link_step_size ) {


    const larcv::ImageMeta& meta = img_v.front().meta();
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick]
    const int nplanes = img_v.size();

    std::vector< std::vector<double> > output_path;

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

	// Give the corresponding values as floats to eliminate the error.
	Float_t xyz_float[3];
	std::vector<float> pt_float(3,0.0);
	
        for (int i=0; i<3; i++) {
          xyz[i] = step0[i] + stepsize*istep*dir3d[i];
          pt[i] = xyz[i];

	  xyz_float[i] = step0[i] + stepsize*istep*dir3d[i];
          pt_float[i] = xyz_float[i];
        }

	// Ensure that the 3 coordinates give a value inside the image before appending this value to the track.
	// Convert the 3D position to wire coordinates (I have to use xyz to do this).
	std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( pt_float, img_v.front().meta(), img_v.size() );

	// Continue if the imgcoords are not inside the image.
	// If any dimension is outside the image, then please continue.

	bool point_gives_2D_pixel_out_of_range = false;
	
	if (imgcoords[0] < 0 || imgcoords[0] >= img_v.at(0).meta().rows())
	  point_gives_2D_pixel_out_of_range = true;

	for (size_t p = 0; p < img_v.size(); p++) {

	  if (imgcoords[p+1] < 0 || imgcoords[p+1] >= img_v.at(0).meta().cols())
	    point_gives_2D_pixel_out_of_range = true;

	}

	if (!point_gives_2D_pixel_out_of_range)
	  output_path.emplace_back( std::move(pt) );

      } //end of steps

    }//end of node loop

    // fill track3d data
    BMTrackCluster3D track3d( start_pt, end_pt, output_path );

    return track3d;
  }

}
