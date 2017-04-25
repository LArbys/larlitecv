#include "UnipolarHackAlgo.h"

#include <cmath>

namespace larlitecv {
  
  std::vector< larcv::Image2D > UnipolarHackAlgo::hackUnipolarArtifact( const std::vector<larcv::Image2D>& img_v, const std::vector<int>& applytoplane,
									const std::vector<float>& neg_threshold ) {
    std::vector< larcv::Image2D > hacked_v;
    for ( size_t p=0; p<img_v.size(); p++) {
      larcv::Image2D hacked( img_v[p] );
      if ( applytoplane[p]>0 ) {
	for (int c=0; c<(int)hacked.meta().cols(); c++) {
	  std::vector<UnipolarROI_t> candidates;
	  scanForPulse( c, hacked, neg_threshold[p], candidates );
	  for ( auto const& candidate : candidates ) {
	    // hack image
	    //std::cout << "fix candidate [" << candidate.end_row << "," << candidate.start_row << " plat_mean=" << candidate.plat_val << "]" << std::endl;
	    // remove plateau
	    for ( int r=candidate.end_row; r<=candidate.start_row; r++ ) {
	      if ( r<0 || r>=hacked.meta().rows() )
		continue;
	      float val = hacked.pixel(r,c);
	      val -= 1.3*candidate.plat_val;
	      //val *= 5.0;
	      if ( val<neg_threshold[p] ) {
		val = 0; // rectify
	      }
	      hacked.set_pixel( r, c, val );
	    }
	    // add hit
	    for (int r=candidate.pos_peak-1; r<=candidate.pos_peak+1; r++){
	      if ( r<0 || r>=hacked.meta().rows() )
		continue;
	      hacked.set_pixel(r,c,80.0);
	    }
	  }
	}
      }
      hacked_v.emplace_back( std::move(hacked) );
    }

    return hacked_v;
  }
  
  void UnipolarHackAlgo::scanForPulse( const int col, const larcv::Image2D& img, const float neg_threshold, std::vector<UnipolarROI_t>& rois ) {
    
    bool inpulse = false;
    bool inplat = false;
    UnipolarROI_t candidate;
    std::deque<float> buff;
    
    for (int r=(int)img.meta().rows()-1; r>=0; r--) {

      float val = img.pixel(r,col);

      if (!inpulse) {
	// while out of peak, look for neg_threshold
	if ( val<neg_threshold ) {
	  // load new candidate
	  inpulse = true;
	  inplat = false;
	  candidate.reset();
	  candidate.neg_peak = r;
	  candidate.pos_peak = r;	  
	  candidate.neg_peak_val = val;
	  candidate.pos_peak_val = val;	  
	  candidate.start_row = r;
	  //std::cout << "  found below-thresh pixel: val=" << val << ". start roi and plateau search."<< std::endl;
	}
      }
      else {
	// in peak
	if ( val>candidate.pos_peak_val ) {
	  candidate.pos_peak_val = val;
	  candidate.pos_peak = r;
	}
	if ( val<candidate.neg_peak_val ) {
	  candidate.neg_peak_val = val;
	  candidate.neg_peak = r;
	}
	
	// while in ROI, we ...
	// 0) update the deq
	buff.push_back( val );
	
	// 1) calculate mean and RMS
	if ( buff.size()>=5 ) {
	  float sum = 0;
	  float sum2 = 0;
	  for (auto const& v : buff) {
	    sum += v;
	    sum2 += v*v;
	  }
	  float mean = sum/buff.size();
	  float mean2 = sum2/buff.size();
	  float var  = sqrt(mean2 - mean*mean);
	  //std::cout << "  add to buff: size=" << buff.size() << " mean=" << mean << " var=" << var << std::endl;
	  //if ( var/mean < 0.1 && mean>4.0 ) {
	  if ( var < 3 && mean>8.0 ) {	  
	    // hit plateau
	    if ( candidate.plat_row<0 ) {
	      // start plateau
	      candidate.plat_row = r;
	      inplat = true;
	      //std::cout << "  in plateau" << std::endl;
	    }
	    if ( inplat ) {
	      candidate.nplat++;
	      candidate.plat_vals.push_back(mean);
	    }
	  }
	  buff.pop_front();
	}//end of if buff.size()>=5

	// see if we've crashed
	if ( inplat && val<2.0 ) {
	  inpulse = false;
	  inplat = false;
	  candidate.end_row = r-1;
	  candidate.plat_val = 0.;
	  for ( auto const& means : candidate.plat_vals )
	    candidate.plat_val += means;
	  candidate.plat_val /= (float)candidate.plat_vals.size();
	  // std::cout << "end of plateau and candidate. define candidate:" << std::endl;
	  // std::cout << "  candidate start tick=" << img.meta().pos_y(candidate.start_row) << std::endl;
	  // std::cout << "  candidate neg peak=" << img.meta().pos_y(candidate.neg_peak) << " val=" << candidate.neg_peak_val << std::endl;
	  // std::cout << "  candidate pos peak=" << img.meta().pos_y(candidate.pos_peak) << " val=" << candidate.pos_peak_val << std::endl;
	  // std::cout << "  candidate plat start=" << img.meta().pos_y( candidate.plat_row ) << std::endl;
	  // std::cout << "  candidate plat mean=" << candidate.plat_val << std::endl;	  
	  // std::cout << "  candidate end tick=" << img.meta().pos_y( candidate.end_row ) << std::endl;	  
	  rois.push_back( candidate );
	}
	  
	      
      }//end of else if in pulse
      
      
    }// loop over rows
  }

}
