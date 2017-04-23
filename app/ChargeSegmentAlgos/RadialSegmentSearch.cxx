#include "RadialSegmentSearch.h"
#include <cstring>
#include <cmath>

#include "DataFormat/ImageMeta.h"

namespace larlitecv {

  std::vector< larcv::Pixel2D > RadialSegmentSearch::findIntersectingChargeClusters( const larcv::Image2D& img, const larcv::Image2D& badch,
										     int row, int col, int radius, float threshold, int arc_divs ) {

    const larcv::ImageMeta& meta = img.meta();
    // first we build an order vector of pixels
    // they need to be unique, so we make an array to mark if we've used a pixel
    
    // we assume (col,row) is still in the image.
    if ( col<0 || col>=(int)meta.cols() || row<0 || row>=(int)meta.rows() || radius<0) {
      std::vector< larcv::Pixel2D > empty;
      return empty;
    }

    // first define a bounding box
    std::vector<int> lowleft(2);
    lowleft[0] = col-radius;
    lowleft[1] = row-radius;
    std::vector<int> upright(2);
    upright[0] = col+radius+1;
    upright[1] = row+radius+1;

    int max[2] = {0};
    max[0] = upright[0]-lowleft[0]+1;
    max[1] = upright[1]-lowleft[1]+1;
    int* markerbox = new int[ max[0]*max[1] ];
    memset(markerbox, 0, sizeof(int));

    std::vector< larcv::Pixel2D > pixelring;
    
    for (int i=0; i<arc_divs; i++) {
      float theta = i*(2*3.14159/float(arc_divs));
      int x = col + radius*cos(theta);
      int y = row + radius*sin(theta);
      if ( x<0 || x>=(int)meta.cols() || y<0 || y>=(int)meta.rows() )
	continue;
      int nx = x-col + radius;
      int ny = y-row + radius;
      if ( markerbox[ max[0]*ny + nx ]==0 ) {
	larcv::Pixel2D pix( x, y );
	if ( badch.pixel(y,x) )
	  pix.Intensity( threshold+1 );
	else
	  pix.Intensity( img.pixel(y,x) );
	pixelring.emplace_back( std::move(pix) );
	markerbox[ max[0]*ny + nx ] = 1;
      }
    }

    delete markerbox;
    
    // loop through and look for pixel spikes
    int idx_start = -1; // pixelering rings

    // first we look for a region, not above threshold to start
    int nbelow_thresh = 0;
    for (size_t i=0; i<pixelring.size(); i++) {
      if ( pixelring[i].Intensity()<threshold )
	nbelow_thresh++;
      else
	nbelow_thresh = 0;
      if ( nbelow_thresh>=3 ) {
	idx_start = (int)i - 3;
	break;
      }
    }

    if ( idx_start<0 ) {
      // everything is above threshold...
      std::vector< larcv::Pixel2D > empty;
      return empty;      
    }

    // now we loop around the ring, finding hits in a very simple manner
    int idx = idx_start+1;

    struct Hit_t {
      int start_idx;
      int end_idx;
      int max_idx;
      float maxval;
      void reset() {
	start_idx = -1;
	end_idx = -1;
	max_idx = -1;
	maxval = -1;
      };
      void update_max( int idx, float val ) {
	if ( val > maxval ) {
	  max_idx = idx;
	  maxval = val;
	}
      };
    };
    std::vector<Hit_t> hitlist;
    
    bool inhit = false;
    while ( idx!=idx_start ) {

      if ( !inhit ) {
	if ( pixelring[idx].Intensity()>=threshold ) {
	  // start a hit
	  Hit_t hit;
	  hit.reset();
	  hit.start_idx = idx;
	  hit.max_idx = idx;
	  hit.maxval = pixelring[idx].Intensity();
	  hitlist.emplace_back( std::move(hit) );
	  inhit = true;
	}
	else {
	  // do nothing
	}
      }
      else if (inhit) {
	if ( pixelring[idx].Intensity()<threshold ) {
	  // end the hit
	  hitlist.back().end_idx = idx-1;
	  inhit = false;
	}
	else {
	  // update the max hit
	  hitlist.back().update_max( idx, pixelring[idx].Intensity() );
	}
      }
      
      idx++;
      // if at end of pixelring vector, loop back around
      if ( idx>=(int)pixelring.size() ) 
	idx = 0;
    }

    delete markerbox;

    // make output pixel2d list. provide location of max of hits along the ring

    std::vector<larcv::Pixel2D> output;
    for ( auto const& hit : hitlist ) {
      larcv::Pixel2D pix( pixelring[ hit.max_idx ] );
      output.emplace_back( std::move(pix) );
    }
    return output;
  }

}
