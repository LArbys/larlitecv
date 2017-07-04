#include "PathPixelChargeMethods.h"

#include <cmath>

#include "DataFormat/Image2D.h"

#include "TaggerTypes/Path2Pixels.h"

namespace larlitecv {

  std::vector< PixelQPt > getPixelQPts( const std::vector< std::vector<double> >& path,
					const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
					const int pixel_radius, const float threshold, const float max_stepsize ) {
    std::vector< PixelQPt > pathq;
    int npts = path.size();
    for ( int ipt=1; ipt<npts; ipt++) {
      const std::vector<double>& curr_pos = path[ipt];
      const std::vector<double>& last_pos = path[ipt-1];

      double norm = 0.;
      std::vector<double> dir(3,0);
      for (int i=0; i<3; i++) {
	dir[i] = curr_pos[i]-last_pos[i];
	norm += dir[i]*dir[i];
      }
      norm = sqrt(norm);
      for (int i=0; i<3; i++)
	dir[i] /= norm;

      int nsubsteps = norm/max_stepsize;
      if ( nsubsteps==0 || (norm - nsubsteps*max_stepsize) > 1.0e-4 )
	nsubsteps++;
      double stepsize = norm/float(nsubsteps);
      if ( ipt==npts-1 )
	nsubsteps++; // for last step, include last substep pos
      for (int iss=0; iss<nsubsteps; iss++) {
	std::vector<double> substep(3,0);
	for (int i=0; i<3; i++)
	  substep[i] = last_pos[i] + dir[i]*(float(iss)*stepsize);
	PixelQPt qpt( substep, img_v, badch_v, pixel_radius, threshold );
	pathq.emplace_back( std::move(qpt) );
      }//end of substep loop
    }//end of step loop

    return pathq;
  }

  std::vector<double> getTrackTotalPixelCharge( const std::vector< std::vector<double> >& path,
						const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
						const int pixel_radius, const float threshold, const float max_stepsize ) {
    std::vector<float> threshold_v( img_v.size(), threshold );
    std::vector<int> neighsize_v( img_v.size(), pixel_radius );
    std::vector<larcv::Pixel2DCluster> track_pixs = getTrackPixelsFromImages( path, img_v, badch_v, threshold_v, neighsize_v, max_stepsize );

    std::vector<double> totalpathq(img_v.size(),0);    
    for ( size_t p=0; p<img_v.size(); p++) {
      for (auto const& pix : track_pixs[p] ) {
	if ( badch_v[p].pixel( pix.Y(), pix.X() )==0 && img_v[p].pixel(pix.Y(),pix.X())>threshold )
	  totalpathq[p] += img_v[p].pixel(pix.Y(),pix.X());
      }
    }
    
    return totalpathq;
  }
  
}
