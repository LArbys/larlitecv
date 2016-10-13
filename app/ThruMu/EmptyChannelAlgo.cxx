#include "EmptyChannelAlgo.h"

// larcv
#include "DataFormat/ImageMeta.h"

namespace larlitecv {

  larcv::Image2D EmptyChannelAlgo::labelEmptyChannels( float threshold, const larcv::Image2D& tpcimg ) {

    // for data I imagine this will need to be more elaborate
    
    const larcv::ImageMeta& meta = tpcimg.meta();
    larcv::Image2D emptyimg( meta );
    emptyimg.paint(0.0);

    for (int col=0; col<meta.cols(); col++) {
      bool isempty = true;
      for ( int row=0; row<meta.rows(); row++) {
	if ( tpcimg.pixel(row,col)>threshold ) {
	  isempty = false;
	  break;
	}
      }

      if ( isempty )
	emptyimg.paint_col( col, threshold );
    }
    
    return emptyimg;
  }

}
