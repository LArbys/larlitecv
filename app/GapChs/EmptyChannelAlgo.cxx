#include "EmptyChannelAlgo.h"
#include <assert.h>

// larcv
#include "DataFormat/ImageMeta.h"
#include "DataFormat/ChStatus.h"

namespace larlitecv {

  larcv::Image2D EmptyChannelAlgo::labelEmptyChannels( float threshold, const larcv::Image2D& tpcimg, const float max_value ) {

    // for data I imagine this will need to be more elaborate
    
    const larcv::ImageMeta& meta = tpcimg.meta();
    larcv::Image2D emptyimg( meta );
    emptyimg.paint(0.0);

    for (int col=0; col<meta.cols(); col++) {
      bool isempty = true;
      for ( int row=0; row<meta.rows(); row++) {
        float val = tpcimg.pixel(row,col);
        if ( val>threshold && (max_value<0 || val<max_value) ) {
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


  std::vector<larcv::Image2D> EmptyChannelAlgo::makeBadChImage( int minstatus, int nplanes, int start_tick, int nticks, int nchannels, 
                                                                int time_downsample_factor, int wire_downsample_factor,
                                                                const larcv::EventChStatus& ev_status ) {
    std::vector<larcv::Image2D> badchs;

    for (int p=0; p<nplanes; p++) {

      larcv::ImageMeta meta( nchannels, nticks, nticks, nchannels, 0.0, start_tick+nticks, p );
      larcv::Image2D img(meta);
      img.paint(0.0);

      auto it_chstatus = ev_status.ChStatusMap().find(p);
      if ( it_chstatus!=ev_status.ChStatusMap().end() ) {
      
        const larcv::ChStatus& chs = ev_status.Status(p);
        const std::vector<short>& status_v = chs.as_vector();
        for ( int ch=0; ch<status_v.size(); ch++) {
          int status = status_v.at(ch);
          if ( ch<nchannels && status<minstatus  )
            img.paint_col( ch, 255.0 );
        }
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

  std::vector<larcv::Image2D> EmptyChannelAlgo::findMissingBadChs( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimgs_v,
    const float empty_ch_threshold, const int max_empty_gap, const float empty_ch_max ) {

    std::cout << "EmptyChannelAlgo::findMissingBadChs" << std::endl;

    std::vector<larcv::Image2D> output;

    for (size_t p=0; p<img_v.size(); p++) {
      const larcv::Image2D& img     = img_v.at(p);
      const larcv::Image2D& badch   = badchimgs_v.at(p);
      const larcv::Image2D& emptych = labelEmptyChannels( empty_ch_threshold, img, empty_ch_max );
      int cols = img.meta().cols(); 
      bool ingap = false;
      int gapsize = 0;
      int gapstart = 0;
      std::vector< std::pair<int,int> > gaplist;
      for ( int c=0; c<cols; c++) {

        if ( !ingap ) {
          // not in a gap
          if ( emptych.pixel(0,c)>0 && badch.pixel(0,c)==0 ) { // we don't want to capture the same info. as the bad channels
            // start a gap
            gapstart = c;
            //if ( gapstart<0 ) gapstart = 0;
            ingap = true;
          }
          else {
            // nothing to do I guess
          }
        }
        else {
          // in a gap
          if ( emptych.pixel(0,c)==0 ) {
            // end a gap
            int start = gapstart;
            int end   = c-1;
            //std::cout << "found gap [" << start << "," << end << "]" << std::endl;
            gaplist.push_back( std::make_pair<int,int>(std::move(start),std::move(end)) );
            gapsize = 0;
            ingap = false;
          }
          else {
            // continue the gap
            gapsize++;
          }
        }
      }//loop over columns

      // ideally we should do something smart, determining gaps based on if we can tell a track is crossing it.
      // for now, we set a gap limit
      larcv::Image2D gapimg(img_v.at(p).meta());
      gapimg.paint(0);

      for ( auto& gap : gaplist ) {
        int gapsize = gap.second-gap.first+1;
        if ( gapsize>max_empty_gap ) continue;
        for (int c=gap.first; c<=gap.second; c++)
          gapimg.paint_col( c, 255.0 );
      }

      output.emplace_back( std::move(gapimg) );

    }//loop over planes

    return output;

  }
  
}
