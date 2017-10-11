#include "Path2Pixels.h"

#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

namespace larlitecv {

  std::vector<larcv::Pixel2DCluster> getTrackPixelsFromImages( const std::vector< std::vector<double> >& path3d,
							       const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
							       const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
							       const float stepsize ) {
    
    if ( thresholds.size()!=imgs.size() || neighborhood_size.size()!=imgs.size() ) {
      throw std::runtime_error("BMTrackCluster3D::getTrackPixelsFromImages[error] number of thresholds or neighborhood size does not match number of images.");
    }
    
    // get unique pixels
    const int nplanes=imgs.size();
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    std::vector< std::set< std::vector<int> > > planepixs(nplanes);
    
    for ( int istep=0; istep<(int)path3d.size()-1; istep++ ) {
      const std::vector<double>& point      = path3d.at(istep);
      const std::vector<double>& next_point = path3d.at(istep+1);
      float dist=0;
      std::vector<double> dir(3);
      for (int i=0; i<3; i++) {
        dir[i] = (next_point[i]-point[i]);
        dist += dir[i]*dir[i];
      }
      dist = sqrt(dist);
      for ( int i=0; i<3; i++ )
        dir[i] /= dist;
      int nsubsteps = dist/stepsize+1;
      float substepsize = dist/float(nsubsteps);
      for (int isub=0; isub<=nsubsteps; isub++) {
        Double_t xyz[3] = {0};
        for (int i=0; i<3; i++) xyz[i] = point[i] + isub*substepsize*dir[i];
        std::vector<float> wires(nplanes);
        bool inplane = true;
        for (int p=0; p<nplanes; p++) {
          wires[p] = larutil::Geometry::GetME()->WireCoordinate(xyz,p);
          if ( wires[p]<imgs.at(p).meta().min_x() || wires[p]>=imgs.at(p).meta().max_x() )
            inplane = false;
        }
        float tick = xyz[0]/cm_per_tick + 3200.0;
        if ( tick<=imgs.front().meta().min_y() || tick>=imgs.front().meta().max_y() )
          inplane = false;
        if ( !inplane ) continue;

        // store the pixel and pixels around it, if above threshold
        int row = imgs.front().meta().row( tick );
        for (int p=0; p<nplanes; p++) {
          int col = imgs.at(p).meta().col(wires[p]);

         for (int dr=-neighborhood_size[p]; dr<=neighborhood_size[p]; dr++) {
           int r = row+dr;
           if ( r<0 || r>=(int)imgs.at(p).meta().rows() ) continue;
           for (int dc=-neighborhood_size[p]; dc<=neighborhood_size[p]; dc++) {
            int c = col+dc;
            if ( c<0 || c>=(int)imgs.at(p).meta().cols() ) continue;
            if ( imgs.at(p).pixel(r,c)>=thresholds.at(p) || badchimgs.at(p).pixel(r,c)>0 ) {
              std::vector<int> pix(2,0);
              pix[0] = c;
              pix[1] = r;
              planepixs[p].insert( pix );
            }
           }
         }
       }
      }//end of substep loop
    }

    // now make pixel clusters
    std::vector<larcv::Pixel2DCluster> pixels(nplanes);

    for ( int p=0; p<nplanes; p++ ) {
      for ( auto const& pix : planepixs[p] ) {
        larcv::Pixel2D pixel( pix[0], pix[1] );
        pixels.at(p).emplace_back( std::move(pixel) );
      }
      //std::cout << "bmtrackcluster3d: plane " << p << " num pixels=" << pixels.at(p).size() << std::endl;
    }
    
    return pixels;
    
  }

  std::vector<larcv::Pixel2DCluster> getTrackPixelsFromImagesNoBadCh( const std::vector< std::vector<double> >& path3d,
								      const std::vector<larcv::Image2D>& imgs, 
								      const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
								      const float stepsize ) {
    
    if ( thresholds.size()!=imgs.size() || neighborhood_size.size()!=imgs.size() ) {
      throw std::runtime_error("BMTrackCluster3D::getTrackPixelsFromImages[error] number of thresholds or neighborhood size does not match number of images.");
    }
    
    // get unique pixels
    const int nplanes=imgs.size();
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    std::vector< std::set< std::vector<int> > > planepixs(nplanes);
    
    for ( int istep=0; istep<(int)path3d.size()-1; istep++ ) {
      const std::vector<double>& point      = path3d.at(istep);
      const std::vector<double>& next_point = path3d.at(istep+1);
      float dist=0;
      std::vector<double> dir(3);
      for (int i=0; i<3; i++) {
        dir[i] = (next_point[i]-point[i]);
        dist += dir[i]*dir[i];
      }
      dist = sqrt(dist);
      for ( int i=0; i<3; i++ )
        dir[i] /= dist;
      int nsubsteps = dist/stepsize+1;
      float substepsize = dist/float(nsubsteps);
      for (int isub=0; isub<=nsubsteps; isub++) {
        Double_t xyz[3] = {0};
        for (int i=0; i<3; i++) xyz[i] = point[i] + isub*substepsize*dir[i];
        std::vector<float> wires(nplanes);
        bool inplane = true;
        for (int p=0; p<nplanes; p++) {
          wires[p] = larutil::Geometry::GetME()->WireCoordinate(xyz,p);
          if ( wires[p]<imgs.at(p).meta().min_x() || wires[p]>=imgs.at(p).meta().max_x() )
            inplane = false;
        }
        float tick = xyz[0]/cm_per_tick + 3200.0;
        if ( tick<=imgs.front().meta().min_y() || tick>=imgs.front().meta().max_y() )
          inplane = false;
        if ( !inplane ) continue;

        // store the pixel and pixels around it, if above threshold
        int row = imgs.front().meta().row( tick );
        for (int p=0; p<nplanes; p++) {
          int col = imgs.at(p).meta().col(wires[p]);

         for (int dr=-neighborhood_size[p]; dr<=neighborhood_size[p]; dr++) {
           int r = row+dr;
           if ( r<0 || r>=(int)imgs.at(p).meta().rows() ) continue;
           for (int dc=-neighborhood_size[p]; dc<=neighborhood_size[p]; dc++) {
            int c = col+dc;
            if ( c<0 || c>=(int)imgs.at(p).meta().cols() ) continue;
            if ( imgs.at(p).pixel(r,c)>=thresholds.at(p) ) {
              std::vector<int> pix(2,0);
              pix[0] = c;
              pix[1] = r;
              planepixs[p].insert( pix );
            }
           }
         }
	}
      }//end of substep loop
    }
    
    // now make pixel clusters
    std::vector<larcv::Pixel2DCluster> pixels(nplanes);
    
    for ( int p=0; p<nplanes; p++ ) {
      for ( auto const& pix : planepixs[p] ) {
        larcv::Pixel2D pixel( pix[0], pix[1] );
        pixels.at(p).emplace_back( std::move(pixel) );
      }
      //std::cout << "bmtrackcluster3d: plane " << p << " num pixels=" << pixels.at(p).size() << std::endl;
    }
    
    return pixels;
    
  }

  std::vector<larcv::Pixel2DCluster> getTrackPixelsFromImages( const larlite::track& lltrack,
							       const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
							       const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
							       const float stepsize ) {
    std::vector<larcv::Pixel2DCluster> pixclust_v;
    // convert larlite track to vector<vector<float>> path. use above methods.
    return pixclust_v;
  }
  
}
