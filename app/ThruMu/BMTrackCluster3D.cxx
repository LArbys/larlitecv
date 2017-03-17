#include "BMTrackCluster3D.h"
#include <cmath>

#include "TVector.h"

namespace larlitecv {

  BMTrackCluster3D::BMTrackCluster3D() {
    start_index = -1;
    end_index = -1;
  }

  BMTrackCluster3D::~BMTrackCluster3D() {}

  void BMTrackCluster3D::markImageWithTrack( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
      const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size, 
      std::vector<larcv::Image2D>& markedimgs, const float markvalue ) const {

    if ( markedimgs.size()==0 ) {
      for (auto &img : imgs ) {
        larcv::Image2D markedimg( img.meta() );
        markedimg.paint(0.0);
        markedimgs.emplace_back( std::move(markedimg) );
      }
    }
      
    if ( thresholds.size()!=imgs.size() || neighborhood_size.size()!=imgs.size() ) {
    	throw std::runtime_error("BMTrackCluster3D::markImageWithTrack[error] number of thresholds or neighborhood size does not match number of images.");
    }

    const std::vector< larlitecv::BMTrackCluster2D >& tracks2d = (*this).plane_paths;
    
    for ( int p=0; p<(int)imgs.size(); p++ ) {

      const larcv::Image2D& img = imgs.at(p);
      const larcv::ImageMeta& meta = img.meta();
      int plane = (int)meta.plane();

      const larlitecv::BMTrackCluster2D& track = tracks2d.at(plane);

      for (size_t ipixel=1; ipixel<track.pixelpath.size(); ipixel++ ) {
        const larcv::Pixel2D& pixel      = track.pixelpath.at(ipixel);
        const larcv::Pixel2D& past_pixel = track.pixelpath.at(ipixel-1);

        int nsteps_row = fabs((float)pixel.Y()-(float)past_pixel.Y());
        int nsteps_col = fabs((float)pixel.X()-(float)past_pixel.X());
        float step[2] = { (float)pixel.X()-(float)past_pixel.X(), (float)pixel.Y()-(float)past_pixel.Y() };

        // how we define the step size depends on which dimension has the largest number of steps
        int nsteps = ( nsteps_row > nsteps_col ) ? nsteps_row : nsteps_col;

        for (size_t i=0; i<2; i++)
          step[i] /= float(nsteps);

        for (int istep=0; istep<=nsteps; istep++ ) {

          int col = past_pixel.X() + istep*step[0];
          int row = past_pixel.Y() + istep*step[1];

          for ( int dc=-neighborhood_size.at(plane); dc<=neighborhood_size.at(plane); dc++ ) {
            for ( int dr=-neighborhood_size.at(plane); dr<=neighborhood_size.at(plane); dr++ ) {
              int r = row+dr;
              int c = col+dc;
              if ( r<0 || r>=(int)meta.rows() ) continue;
              if ( c<0 || c>=(int)meta.cols() ) continue;
              float val = img.pixel( r, c );
              if ( val> thresholds.at(plane) || badchimgs.at(plane).pixel(r,c)>0 )
                markedimgs.at(p).set_pixel(r,c,markvalue);
            }//end of loop over r-neighborhood
          }//end of loop over c-neighborhood
        }//end of steps
      }//end of pixel list loop

    }//end of img loop
    
  }

  larlite::track BMTrackCluster3D::makeTrack() const {

    larlite::track lltrack;
    int istep = 0;
    for ( auto const& point3d : path3d ) {
      TVector3 vec( point3d[0], point3d[1], point3d[2] );
      lltrack.add_vertex( vec );
      if ( istep+1<(int)path3d.size() ) {
        TVector3 dir( path3d.at(istep+1)[0]-point3d[0], path3d.at(istep+1)[1]-point3d[1], path3d.at(istep+1)[2]-point3d[2] );
        lltrack.add_direction( dir );
      }
      else {
        TVector3 dir( point3d[0]-path3d.at(istep-1)[0], point3d[1]-path3d.at(istep-1)[1], point3d[2]-path3d.at(istep-1)[2] );
        lltrack.add_direction( dir );
      }
    }
    return lltrack;
  }

}
