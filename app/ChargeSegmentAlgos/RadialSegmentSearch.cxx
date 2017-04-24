#include "RadialSegmentSearch.h"
#include <cstring>
#include <cmath>

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

#include "DataFormat/ImageMeta.h"

#include "RadialSegmentSearchTypes.h"

namespace larlitecv {

  std::vector< larlitecv::RadialHit_t > RadialSegmentSearch::findIntersectingChargeClusters( const larcv::Image2D& img, const larcv::Image2D& badch,
											     int row, int col, int radius, float threshold, int arc_divs ) {

    const larcv::ImageMeta& meta = img.meta();
    // first we build an order vector of pixels
    // they need to be unique, so we make an array to mark if we've used a pixel
    
    // we assume (col,row) is still in the image.
    if ( col<0 || col>=(int)meta.cols() || row<0 || row>=(int)meta.rows() || radius<0) {
      std::vector< larlitecv::RadialHit_t > empty;
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
    memset(markerbox, 0, sizeof(int)*max[0]*max[1] );

    std::cout << "made pixel box: upright(" << upright[0] << "," << upright[1] << ") to lowerleft(" << lowleft[0] << "," << lowleft[1] << ")" << std::endl;
    
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
    std::cout << "pixels in ring: " << pixelring.size() << std::endl;

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
	idx_start = (int)i-3+1;
	break;
      }
    }

    std::cout << "starting ring index: " << idx_start << std::endl;
    if ( idx_start<0 ) {
      // everything is above threshold...
      std::vector< larlitecv::RadialHit_t > empty;
      return empty;      
    }

    // now we loop around the ring, finding hits in a very simple manner
    int idx = idx_start+1;

    std::vector< larlitecv::RadialHit_t> hitlist;
    
    bool inhit = false;
    bool loopedaround = false;
    while ( idx!=idx_start ) {
      if ( !inhit ) {
	if ( pixelring[idx].Intensity()>=threshold ) {
	  // start a hit
	  RadialHit_t hit;
	  hit.reset();
	  hit.start_idx = idx;
	  hit.max_idx = idx;
	  hit.maxval = pixelring[idx].Intensity();
	  hit.pixlist.push_back( pixelring[idx] );
	  hitlist.emplace_back( std::move(hit) );
	  inhit = true;
	  std::cout << "start hit at " << idx << std::endl;
	}
	else {
	  // do nothing
	}
      }
      else if (inhit) {
	if ( pixelring[idx].Intensity()<threshold ) {
	  // end the hit
	  if ( !loopedaround )
	    hitlist.back().set_end( idx-1 );
	  else
	    hitlist.back().set_end( idx-1+pixelring.size() );
	  
	  inhit = false;
	  std::cout << "end hit at " << idx << std::endl;
	}
	else {
	  // update the max hit
	  if ( !loopedaround )
	    hitlist.back().update_max( idx, pixelring[idx].Intensity() );
	  else
	    hitlist.back().update_max( idx+pixelring.size(), pixelring[idx].Intensity() );
	  hitlist.back().add_pixel( pixelring[idx] );
	}
      }
      idx++;
      // if at end of pixelring vector, loop back around
      if ( idx>=(int)pixelring.size() )  {
	idx = 0;
	loopedaround = true;
      }
    }

    return hitlist;
  }

  std::vector< larlitecv::RadialHit_t > RadialSegmentSearch::findIntersectingChargeClustersX( const larcv::Image2D& img, const larcv::Image2D& badch,
											      const std::vector<float>& pos3d, float radius, float threshold ) {

    const larcv::ImageMeta& meta = img.meta();
    // first we build an order vector of pixels
    // they need to be unique, so we make an array to mark if we've used a pixel

    std::vector< std::vector<Double_t> > yzbox;
    for (int z=0; z<2; z++) {
      for (int y=0; y<2; y++) {
	std::vector<Double_t> yz(3,0.0);
	yz[1] = pos3d[1] + (2*y-1)*radius;
	yz[2] = pos3d[2] + (2*z-1)*radius;
	yzbox.push_back(yz);
      }
    }

    std::vector<int> colbounds(2,-1);
    for ( size_t i=0; i<yzbox.size(); i++ ) {
      float wire = larutil::Geometry::GetME()->WireCoordinate( yzbox[i], (int)meta.plane() );
      wire = ( wire<0 ) ? 0 : wire;
      wire = (wire>=meta.max_x()) ? meta.max_x()-1 : wire;
      wire = (int)meta.col(wire);
      
      if ( colbounds[0]<0 || colbounds[0]>wire )
	colbounds[0] = wire;
      if ( colbounds[1]<0 || colbounds[1]<wire )
	colbounds[1] = wire;
    }

    float tick = pos3d[0]/(larutil::LArProperties::GetME()->DriftVelocity()*0.5)+3200.0;
    int row = meta.row( tick );
    int rowradius = radius/(larutil::LArProperties::GetME()->DriftVelocity()*0.5)/meta.pixel_height();
    
    // first define a bounding box
    std::vector<int> lowleft(2);
    lowleft[0] = colbounds[0];
    lowleft[1] = row-rowradius;
    lowleft[1] = ( lowleft[1]<0 ) ? 0 : lowleft[1];
    std::vector<int> upright(2);
    upright[0] = colbounds[1]+1;
    upright[1] = row+rowradius+1;
    upright[1] = ( upright[1]>=(int)meta.rows() ) ? (int)meta.rows()-1 : upright[1];

    std::cout << "made pixel box: upright(" << upright[0] << "," << upright[1] << ") to lowerleft(" << lowleft[0] << "," << lowleft[1] << ")" << std::endl;
    
    std::vector< larcv::Pixel2D > pixelring;
    // make a box of pixels
    for (int c=colbounds[0]; c<=colbounds[1]; c++) {
      larcv::Pixel2D pix( c, upright[1] );
      pix.Intensity( img.pixel(upright[1],c) );
      pixelring.push_back( pix );
    }
    for (int r=upright[1]; r>=lowleft[1]; r--) {
      larcv::Pixel2D pix( colbounds[1], r );
      pix.Intensity( img.pixel(r,colbounds[1]) );
      pixelring.push_back( pix );
    }
    for (int c=colbounds[1]; c>=colbounds[0]; c--) {
      larcv::Pixel2D pix( c, lowleft[1] );
      pix.Intensity( img.pixel(lowleft[1],c) );
      pixelring.push_back( pix );
    }
    for (int r=lowleft[1]; r<=upright[1]; r++) {
      larcv::Pixel2D pix( colbounds[0], r );
      pix.Intensity( img.pixel(r,colbounds[0]) );
      pixelring.push_back( pix );
    }
    
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
	idx_start = (int)i-3+1;
	break;
      }
    }

    std::cout << "starting ring index: " << idx_start << std::endl;
    if ( idx_start<0 ) {
      // everything is above threshold...
      std::vector< larlitecv::RadialHit_t > empty;
      return empty;      
    }

    // now we loop around the ring, finding hits in a very simple manner
    int idx = idx_start+1;

    std::vector< larlitecv::RadialHit_t> hitlist;
    
    bool inhit = false;
    bool loopedaround = false;
    while ( idx!=idx_start ) {
      if ( !inhit ) {
	if ( pixelring[idx].Intensity()>=threshold ) {
	  // start a hit
	  RadialHit_t hit;
	  hit.reset();
	  hit.start_idx = idx;
	  hit.max_idx = idx;
	  hit.maxval = pixelring[idx].Intensity();
	  hit.pixlist.push_back( pixelring[idx] );
	  hitlist.emplace_back( std::move(hit) );
	  inhit = true;
	  std::cout << "start hit at " << idx << std::endl;
	}
	else {
	  // do nothing
	}
      }
      else if (inhit) {
	if ( pixelring[idx].Intensity()<threshold ) {
	  // end the hit
	  if ( !loopedaround )
	    hitlist.back().set_end( idx-1 );
	  else
	    hitlist.back().set_end( idx-1+pixelring.size() );
	  
	  inhit = false;
	  std::cout << "end hit at " << idx << std::endl;
	}
	else {
	  // update the max hit
	  if ( !loopedaround )
	    hitlist.back().update_max( idx, pixelring[idx].Intensity() );
	  else
	    hitlist.back().update_max( idx+pixelring.size(), pixelring[idx].Intensity() );
	  hitlist.back().add_pixel( pixelring[idx] );
	}
      }
      idx++;
      // if at end of pixelring vector, loop back around
      if ( idx>=(int)pixelring.size() )  {
	idx = 0;
	loopedaround = true;
      }
    }

    return hitlist;
  }
  
}
