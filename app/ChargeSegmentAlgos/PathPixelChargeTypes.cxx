#include "PathPixelChargeTypes.h"

#include <sstream>

#include "UBWireTool/UBWireTool.h"


namespace larlitecv {

  PixelQPt::PixelQPt( const std::vector<double>& pos3d, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const int pixel_radius, const float threshold ) {
    // we are given a 3d position and images.
    // we project the point back into the images, sum the charge in a box around the center pixel with pixel_radius
    // we also store the meta info as well. (not now, don't really need it)

    // check variables
    if ( pos3d.size()!=3 ) {
      std::stringstream msg;
      msg << __FILE__ << ":" << __LINE__ << " pos3d vector is not size 3";
      throw std::runtime_error(msg.str());
    }
    if ( pixel_radius<0 ) {
      std::stringstream msg;
      msg << __FILE__ << ":" << __LINE__ << " pixel_radius is negative.";
      throw std::runtime_error(msg.str());
    }
    
    std::vector<float> fpos( pos3d.size(), 0); // types suck
    m_pos3d.resize(3,0);
    for (int i=0; i<3; i++) {
      fpos[i] = pos3d[i];
      m_pos3d[i] = pos3d[i];
    }
    std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( fpos, img_v.front().meta(), img_v.size() );
    m_imagerow = imgcoords[0];
    m_imagecols.resize( img_v.size(), 0 );
    for (size_t p=0; p<img_v.size(); p++)
      m_imagecols[p] = imgcoords[p+1];

    m_planeq.resize( img_v.size(), 0 );
    m_goodpixels_per_plane.resize( img_v.size(), 0 );
    m_badchpixels_per_plane.resize( img_v.size(), 0 );
    for (int r=-pixel_radius; r<=pixel_radius; r++) {
      int row = m_imagerow+r;
      if ( row<0 || row>=img_v.front().meta().rows() )
	continue;
      for (size_t p=0; p<img_v.size(); p++) {
	for (int c=-pixel_radius; c<=pixel_radius; c++) {
	  int col = m_imagecols[p] + c;	  
	  if ( col<0 || col>=img_v.front().meta().cols() )
	    continue;

	  // get charge (or badch)
	  float q = img_v[p].pixel(row,col);
	  if ( q>threshold ) {
	    m_planeq[p] += q;
	    m_goodpixels_per_plane[p]++;
	  }
	  else if ( badch_v[p].pixel(row,col)>0 ) {
	    m_badchpixels_per_plane[p]++;
	  }
	}//end of col loop
      }// end of plane loop
    }//end of r loop
  }//end of constructor
  
}
