#include "EmptyChannelAlgo.h"
#include <assert.h>

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

  std::vector<larcv::Image2D> EmptyChannelAlgo::makeBadChImage( int minstatus, int nplanes, int start_tick, int nticks, int nchannels, 
								int time_downsample_factor, int wire_downsample_factor,
								const larlite::event_chstatus& ev_status ) {
    std::vector<larcv::Image2D> badchs;
    for (int p=0; p<nplanes; p++) {

      if ( p>=ev_status.size() ) {
	std::cout << "ch status not available for plane " << p << std::endl;
	assert(false);
      }

      larcv::ImageMeta meta( nchannels, nticks, nticks, nchannels, 0.0, start_tick+nticks, p );
      larcv::Image2D img(meta);
      img.paint(0.0);
      
      const larlite::chstatus& chs = ev_status.at(p);
      const std::vector<short>& status_v = chs.status();
      for ( int ch=0; ch<status_v.size(); ch++) {
	int status = status_v.at(ch);
	if ( ch<nchannels && status<minstatus  )
	  img.paint_col( ch, 255.0 );
      }
      
      if (time_downsample_factor!=1 || wire_downsample_factor!=1 ) {
	int new_rows = meta.rows()/time_downsample_factor;
	int new_cols = meta.cols()/wire_downsample_factor;
	img.compress( new_rows, new_cols );
      }

      badchs.emplace_back( img );
    }
    return badchs;
  }
}
