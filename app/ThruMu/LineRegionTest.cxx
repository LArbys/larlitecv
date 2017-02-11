#include "LineRegionTest.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace larlitecv {

  bool LineRegionTest::test( const BoundarySpacePoint& start_v, const BoundarySpacePoint& end_v,
                             const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
                             std::vector< BMTrackCluster2D >* opt_track2d ) {
    // simply make a line through the end points
    // check in the neighborhood of that line if we have a line of charge going from one to the other

    bool pass_plane[3] = {false, false, false};
    for (int p=0; p<3; p++) {
      // get direction between start and end
      float dir[2] = { (float)(end_v.at(p).col-start_v.at(p).col), (float)(end_v.at(p).row-start_v.at(p).row) };
      int gap[2] =   { (int)(end_v.at(p).col-start_v.at(p).col), (int)(end_v.at(p).row-start_v.at(p).row) };
      // determine which image dimension we loop over and which one we define a region over
      int loopidx = 0; // loop over wire (assuming wire direction component bigger)
      int regionidx = 1; // check regions over time
      if ( fabs(dir[1])>fabs(dir[0]) ) {
        // or vice versa if time direction component bigger
        loopidx = 1; // loop around col axis
        regionidx = 0; // check regions over wires
      }

      // normalize direction
      float norm = 0;
      for (int i=0; i<2; i++)
        norm += dir[i]*dir[i];
      norm = sqrt(norm);
      for (int i=0; i<2; i++)
        dir[i] /= norm;

      float stepsize[2] = {0,0};
      stepsize[loopidx] = dir[loopidx]/fabs(dir[loopidx]);
      stepsize[regionidx] = dir[regionidx]/fabs(dir[loopidx]);
      // note: because we check tha loopidx dir bigger, we should not have a case where we divide by zero

      int nsteps = abs(gap[loopidx]);
      int npass = 0;
      int ncheck = 0;
      float last_pos[2] = { 0., 0.}; // (c,r)
      if ( opt_track2d!=NULL ) {
        larcv::Pixel2D node( start_v.at(p).getcol(), start_v.at(p).getrow() );
        node.Intensity(0);
        opt_track2d->at(p).pixelpath.emplace_back( std::move(node) );
      }
      for (int istep=1; istep<nsteps; istep++) {
        int pos[2] =  { (int)(start_v.at(p).getcol() + istep*stepsize[0]), (int)(start_v.at(p).getrow() + istep*stepsize[1]) };
        if ( verbose_debug )
          std::cout << "LRT step=" << istep << " pos=(" << pos[0] << "," << pos[1] << ")";
        if ( pos[0]<0 ||  pos[0]>=(int)img_v.at(p).meta().cols() || pos[1]<0 || pos[1]>=(int)img_v.at(p).meta().rows() )  {
          if ( verbose_debug ) std::cout << " out of image." << std::endl;
          continue;
        }

        ncheck++;
        bool found_charge = false;
        bool badch = false;
        int maxidx = 0;
        float maxval = 0;
        if ( regionidx==1 ) {
          // region scan in time dimension
          for (int r=-fRegionWidth;r<=fRegionWidth; r++) {
            if ( pos[1]+r>=0 && pos[1]+r<(int)img_v.at(p).meta().rows()
                 && ( (int)img_v.at(p).pixel(pos[1]+r,pos[0])>fPixelThreshold || badchimg_v.at(p).pixel(pos[1]+r,pos[0])>0 ) ) {
              found_charge = true;
              if ( opt_track2d==NULL )
                break;
              else {
                float val = img_v.at(p).pixel(pos[1]+r,pos[0]);
                if ( val>maxval ) {
                  maxidx = pos[1]+r;
                  maxval = val;
                }
                if ( val<fPixelThreshold && badchimg_v.at(p).pixel(pos[1]+r,pos[0])>0 ) {
                  badch = true;
                }
              }
            }
          }
          if ( opt_track2d!=NULL && !badch ) {
            last_pos[0] = pos[0];
            last_pos[1] = maxidx;
            larcv::Pixel2D node( last_pos[0], last_pos[1] );
            opt_track2d->at(p).pixelpath.emplace_back( std::move(node) );
            node.Intensity( maxidx-pos[1] );
          }
        }//end of if regionidx==1
        else {
          /// region scan in column direction
          for (int c=-fRegionWidth;c<=fRegionWidth; c++) {
            if ( pos[0]+c>=0 && pos[0]+c<(int)img_v.at(p).meta().cols()
                 && ( img_v.at(p).pixel(pos[1],pos[0]+c)>fPixelThreshold ) ) {
              found_charge = true;
              if ( opt_track2d==NULL )
                break;
              else {
                float val = img_v.at(p).pixel(pos[1],pos[0]+c);
                if ( val>maxval ) {
                  maxidx = pos[0]+c;
                  maxval = val;
                }
                if ( val<fPixelThreshold && badchimg_v.at(p).pixel(pos[1],pos[0]+c)>0 ) {
                  badch = true;
                }
              }
            }
          }
          if ( pos[0]>=0 && pos[0]<(int)img_v.at(p).meta().cols() && badchimg_v.at(p).pixel(pos[1],pos[0])>0 ) {
            // badch test (a little different from above)
            float val = img_v.at(p).pixel(pos[1],pos[0]);
            if ( val<fPixelThreshold ) {
              badch = true;
              found_charge = true;
              maxidx = pos[0];
              maxval = val;
            }
          }  
          if ( opt_track2d!=NULL && !badch ) {
            last_pos[0] = maxidx;
            last_pos[1] = pos[1];
            larcv::Pixel2D node( last_pos[0], last_pos[1] );
            opt_track2d->at(p).pixelpath.emplace_back( std::move(node) );
            node.Intensity( maxidx-pos[0] );
          }
        }
        if ( found_charge ) {
          npass++;
          if ( verbose_debug ) {
            if ( !badch )
              std::cout << " found charge. maxidx=" << maxidx << std::endl;
            else
              std::cout << " badch" << std::endl;
          }
        }
        else {
          if ( verbose_debug )
            std::cout << " no charge." << std::endl;
        }
      }//end of loopover steps
      
      if ( opt_track2d!=NULL ) {
        larcv::Pixel2D node( end_v.at(p).getcol(), end_v.at(p).getrow() );
        node.Intensity(0);
        opt_track2d->at(p).pixelpath.emplace_back( std::move(node) );
      }

      last_fractions[p] = float(npass)/float(ncheck);
      if (  last_fractions[p] > fFractionThreshold )
        pass_plane[p] = true;
      
    }//end of loop over planes
    
    int npass = 0;
    for (int p=0; p<3; p++) {
      if ( pass_plane[p] ) npass++;
    }
    
    if ( npass>=2 )
      return true;

    return false;

  }


}
