#include "Segment3DAlgo.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv
#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  std::vector< Segment3D_t > Segment3DAlgo::find3DSegments( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
							    const int row_a, const int row_b, const std::vector<float>& thresholds, const int min_hit_width, const int hit_neighborhood ) {

    // first find a list of "hits" for each row on each plane
    int lowrow = row_a;
    int highrow = row_b;
    if ( lowrow > highrow ) {
      lowrow = row_b;
      highrow = row_a;
    }

    std::vector< std::vector<Segment2D_t> > plane_segments2d;
    for (size_t p=0; p<img_v.size(); p++) {
      // get hits
      std::vector<int> hits_low  = findHits( img_v.at(p), lowrow,  thresholds[p], min_hit_width );
      std::vector<int> hits_high = findHits( img_v.at(p), highrow, thresholds[p], min_hit_width );

      // get 2D line segments
      std::vector< Segment2D_t > segments2d = make2DSegments( img_v.at(p), badch_v.at(p), lowrow, hits_low, highrow, hits_high, thresholds[p], hit_neighborhood, 0.9 );
      plane_segments2d.emplace_back( std::move(segments2d) );
    }

    std::vector< Segment3D_t > segments;
    combine2Dinto3D( plane_segments2d, img_v, badch_v, hit_neighborhood, thresholds, 0.9, segments );

    return segments;
  }

  std::vector<int> Segment3DAlgo::findHits( const larcv::Image2D& img, const int row, const float& threshold, const int min_hit_width ) {
    std::vector<int> hits;

    // simple on-off hit finder with min with

    // start off, unless above threshold
    bool onstate = false;
    int hitwidth = 0;
    float maxval = -1.;
    int maxcol = -1;
    if ( img.pixel( row, 0 )>threshold ) {
      onstate = true;
      maxval = img.pixel(row,0);
      maxcol = 0;
      hitwidth++;
    }

    const larcv::ImageMeta& meta = img.meta();

    for (int c=0; c<(int)meta.cols(); c++) {
      float val = img.pixel(row,c);
      if ( onstate ) {
	if ( val>=threshold ) {
	  if ( maxval<val ){
	    maxval = val;
	    maxcol = c;
	  }
	  hitwidth++;
	}
	else if ( val<threshold ) {
	  onstate = false;
	  if ( hitwidth>=min_hit_width ) {
	    // define a new hit!
	    int hitcenter = maxcol;
	    hits.push_back(hitcenter);
	  }
	  hitwidth = 0;
	  maxval = -1.0;
	  maxcol = -1;
	}
      }
      else {
	if ( val<threshold )
	  continue;
	else {
	  // start a hit
	  onstate = true;
	  maxval = val;
	  maxcol = c;
	  hitwidth++;
	}
      }
    }//end of col loop

    return hits;
  }

  std::vector< Segment2D_t > Segment3DAlgo::make2DSegments( const larcv::Image2D& img, const larcv::Image2D& badch, const int lowrow, const std::vector<int>& hits_low,
							    const int highrow, const std::vector<int>& hits_high, const float threshold,
							    const int hit_width, const float frac_good ) {

    //std::cout << "make2DSegments" << std::endl;
    std::vector<Segment2D_t> segments;

    for (int ilow=0; ilow<(int)hits_low.size(); ilow++) {
      for (int ihigh=0; ihigh<(int)hits_high.size(); ihigh++) {

	int low  = hits_low[ilow];
	int high = hits_high[ihigh];

	int num_rows, nrows_w_charge;
	checkSegmentCharge( img, badch, lowrow, low, highrow, high, hit_width, threshold, nrows_w_charge, num_rows );
	float frac = nrows_w_charge/float(num_rows);
	//std::cout << "hit (lo,hi) combo: (" << low << "," << high << ") "
	//	  << "nrows w/ q=" << nrows_w_charge << " nrows=" << num_rows << " frac_good=" << frac << std::endl;

	if ( frac>frac_good ) {
	  Segment2D_t seg2d;
	  seg2d.row_high = highrow;
	  seg2d.row_low  = lowrow;
	  seg2d.col_high = high;
	  seg2d.col_low  = low;
	  segments.emplace_back( std::move(seg2d) );
	}
      }// high loop
    }// low loop

    return segments;
  }

  void Segment3DAlgo::checkSegmentCharge( const larcv::Image2D& img, const larcv::Image2D& badch, const int low_row, const int low_col, const int high_row, const int high_col,
					  const int hit_neighborhood, const float threshold, int& nrows_w_charge, int& num_rows ) {
    //std::cout << "checkSegmentCharge" << std::endl;
    int drow  = abs(high_row-low_row);
    int dcol  = high_col-low_col;
    float slope = float(dcol)/float(drow);
    
    bool hascharge = false;
    nrows_w_charge = 0;
    num_rows = 0;
    for (int r=0; r<=drow; r++) {
      int pixr = low_row + r;
      int pixc = low_col + slope*r;
      if ( pixr<0 || pixr>=(int)img.meta().rows()) continue;
      hascharge = false;
      for (int dr=-hit_neighborhood; dr<=hit_neighborhood; dr++) {
	int rtest = pixr + dr;
	if ( rtest<0 || rtest>=(int)img.meta().rows()) continue;
	for (int dc=-hit_neighborhood; dc<=hit_neighborhood; dc++) {
	  int ctest = pixc+dc;
	  if ( ctest<0 || ctest>=(int)img.meta().cols() ) continue;
	  if ( img.pixel(rtest,ctest)>threshold || badch.pixel(rtest,ctest)>0) {
	    hascharge = true;
	    break;
	  }
	}
      }
      //std::cout << " check (c,r)=(" << pixc << "," << pixr << ") hascharge=" << hascharge << std::endl;
      if ( hascharge )
	nrows_w_charge++;
      num_rows++;
    }//end of row loop
    
  }


  void Segment3DAlgo::combine2Dinto3D( const std::vector< std::vector<Segment2D_t> >& plane_segments2d, const std::vector< larcv::Image2D >& img_v, const std::vector<larcv::Image2D>& badch_v,
				       const int hit_width, const std::vector<float>& thresholds, float good_frac, std::vector<Segment3D_t>& segments ) {

    // We get combinations of Segment2D_t and check their 3 plane consistency.
    // it is assumed that the segments are over the same time interval

    struct SegID {
      size_t plane;
      size_t idx;
    };

    int num_segs = 0;
    for ( auto const& pseg : plane_segments2d )
      num_segs += (int)pseg.size();
    std::vector<SegID> segidx( num_segs );

    int idx=0;
    for (size_t p=0; p<plane_segments2d.size(); p++) {
      for ( size_t i=0; i<plane_segments2d.at(p).size(); i++ ) {
        segidx[idx].plane = p;
        segidx[idx].idx   = i;
        idx++;
      }
    }

    if ( verbosity>1 )
      std::cout << "    Segment3DAlgo: Combine2d into 3d" << std::endl;

    for (int i=0; i<num_segs; i++) {
      for (int j=i+1; j<num_segs; j++) {
        // don't compare 2D segments on the same plane
        if ( segidx[i].plane==segidx[j].plane )
          continue;
        const SegID& idx_a = segidx[i];
        const SegID& idx_b = segidx[j];

        const Segment2D_t& seg_a = plane_segments2d.at(idx_a.plane).at(idx_a.idx);
        const Segment2D_t& seg_b = plane_segments2d.at(idx_b.plane).at(idx_b.idx);

        // what regime are we dealing with
        bool same_high_row = (seg_a.row_high==seg_b.row_high);
        bool same_low_row  = (seg_a.row_low==seg_b.row_low);

        int high_wire_a = img_v.at(idx_a.plane).meta().pos_x( seg_a.col_high );
        int low_wire_a  = img_v.at(idx_a.plane).meta().pos_x( seg_a.col_low );
        int high_wire_b = img_v.at(idx_b.plane).meta().pos_x( seg_b.col_high );
        int low_wire_b  = img_v.at(idx_b.plane).meta().pos_x( seg_b.col_low );
        int row_high = seg_a.row_high;
        int row_low  = seg_a.row_low;

        if ( same_high_row && same_low_row ) {
          //made in the shade
	  if ( verbosity>1 )
	    std::cout << "      same low (" << seg_a.row_low << ") and high (" << seg_a.row_high << ") rows." << std::endl;
        }
        else if ( same_high_row ) {
          // we have to adjust the bottom ends to match at same row
          int target_row = (seg_a.row_low<seg_b.row_low) ? seg_a.row_low : seg_b.row_low;
          float dcol_a = (low_wire_a-high_wire_a);
          float drow_a = (seg_a.row_low-seg_a.row_high);
          float ftarget_wire_a = ( drow_a!=0 ) ? (dcol_a/drow_a)*(target_row-seg_a.row_high) + high_wire_a : low_wire_a;
          float dcol_b = (low_wire_b-high_wire_b);
          float drow_b = (seg_b.row_low-seg_b.row_high);
          float ftarget_wire_b = ( drow_b!=0 ) ? (dcol_b/drow_b)*(target_row-seg_b.row_high) + high_wire_b : low_wire_b;
          low_wire_a = ftarget_wire_a;
          low_wire_b = ftarget_wire_b;
          row_low = target_row;
	  if ( verbosity>1 )  {
	    std::cout << "      same high (" << seg_a.row_high << ") rows. adjusted low= " << row_low
		      << " A(p=" << segidx[i].plane << "): (" << seg_a.row_low << "," << seg_a.col_low << ")->(" << row_low << "," << low_wire_a << ")"
		      << " B(p=" << segidx[j].plane << "): (" << seg_b.row_low << "," << seg_b.col_low << ")->(" << row_low << "," << low_wire_b << ")"
		      << std::endl;
	  }
	  
        }
        else if ( same_low_row ) {
          // we adjust so that top ends match at same row
          int target_row = (seg_a.row_high>seg_b.row_high) ? seg_a.row_high : seg_b.row_high;
          float dcol_a = (high_wire_a-low_wire_a);
          float drow_a = (seg_a.row_high-seg_a.row_low);
          float ftarget_wire_a = ( drow_a!=0 ) ? (dcol_a/drow_a)*fabs(target_row-seg_a.row_low) + low_wire_a : high_wire_a;
          float dcol_b = (high_wire_b-low_wire_b);
          float drow_b = (seg_b.row_high-seg_b.row_low);
          float ftarget_wire_b = ( drow_b!=0 ) ? (dcol_b/drow_b)*fabs(target_row-seg_b.row_low) + low_wire_b : high_wire_b;
          high_wire_a = ftarget_wire_a;
          high_wire_b = ftarget_wire_b;
          row_high = target_row;
	  if ( verbosity>1 )  {
	    std::cout << "      same low (" << seg_a.row_low << ") rows. adjusted high= " << row_high
		      << " low_wire_a=" << high_wire_a << " low_wire_b=" << high_wire_b
		      << std::endl;
	  }
	  
        }
        else {
	  std::cout << "    should not get here." << std::endl;
          continue;
        }

        std::vector<float> poszy_high;
        int crosses_high = 0;
        int otherwire_high = 0;
        int otherplane_high = 0;
        larcv::UBWireTool::getMissingWireAndPlane( (int)idx_a.plane, high_wire_a, (int)idx_b.plane, high_wire_b, otherplane_high, otherwire_high, poszy_high, crosses_high );
        if ( crosses_high==0 ) {
          //std::cout << "    high hit does not cross." << std::endl;
          continue; // these wires don't intersect anywhere
        }

        std::vector<float> poszy_low;
        int crosses_low = 0;
        int otherwire_low = 0;
        int otherplane_low = 0;
        larcv::UBWireTool::getMissingWireAndPlane( (int)idx_a.plane, low_wire_a, (int)idx_b.plane, low_wire_b, otherplane_low, otherwire_low, poszy_low, crosses_low );
        if ( crosses_low==0 ) {
          //std::cout << "    low hit does not cross." << std::endl;
          continue; // these wires don't intersect anywhere
        }

        // ok so both wires cross.
        if ( otherwire_high<0 || otherwire_high>=img_v.at(otherplane_high).meta().max_x() ) {
          //std::cout << "    high is out of bounds" << std::endl;
          continue;
        }
        if ( otherwire_low<0 || otherwire_low>=img_v.at(otherplane_low).meta().max_x() ) {
          //std::cout << "     low is out of bounds" << std::endl;
          continue;
        }
        // let's check the other plane
        int othercol_high = img_v.at(otherplane_high).meta().col( otherwire_high );
        int othercol_low  = img_v.at(otherplane_low).meta().col( otherwire_low );

        int nrows_w_charge = 0;
        int nrows;
        checkSegmentCharge( img_v.at(otherplane_high), badch_v.at(otherplane_high), seg_a.row_low, othercol_low, seg_a.row_high, othercol_high,
			    hit_width, thresholds.at(otherplane_high), nrows_w_charge, nrows );
        // std::cout << "     2d combo(" << i << "," << j << ") "
	// 	  << " rowhi=" << row_high << " (p" << idx_a.plane << "," << high_wire_a << ") (p" << idx_b.plane << "," << high_wire_b << ") "
	// 	  << " rowlo=" << row_low  << " (p" << idx_a.plane << "," << low_wire_a  << ") (p" << idx_b.plane << "," << low_wire_b << ") "
	// 	  << " otherplane=" << otherplane_high << " nrows=" << nrows << " frac=" << float(nrows_w_charge)/nrows << std::endl;
        if ( nrows>0 && float(nrows_w_charge)/nrows > good_frac ) {
          // make segment3d
          Segment3D_t seg3d;
          seg3d.start[1] = poszy_high[1];
          seg3d.start[2] = poszy_high[0];
          seg3d.start[0] = (img_v.at(otherplane_high).meta().pos_y( row_high )-3200.0)*(larutil::LArProperties::GetME()->DriftVelocity()*0.5);
          seg3d.end[1]   = poszy_low[1];
          seg3d.end[2]   = poszy_low[0];
          seg3d.end[0]   = (img_v.at(otherplane_low).meta().pos_y(  row_low )-3200.0)*(larutil::LArProperties::GetME()->DriftVelocity()*0.5);

          // check it's not repeating
          bool isduplicate = false;
          for ( auto const& pastseg : segments ) {
            float dist_hi = 0.;
            float dist_lo = 0;
            for (int i=0; i<3; i++) {
              dist_hi += (seg3d.start[i]-pastseg.start[i])*(seg3d.start[i]-pastseg.start[i]);
              dist_lo += (seg3d.end[i]-pastseg.end[i])*(seg3d.end[i]-pastseg.end[i]);
            }
            dist_hi = sqrt(dist_hi);
            dist_lo = sqrt(dist_lo);
            //std::cout << "     dist_hi=" << dist_hi << " dist_lo=" << dist_lo << std::endl;
            if ( dist_hi<0.8 && dist_lo<0.8 ) {
              isduplicate = true;
              break;
            }
          }

          if ( !isduplicate ) {
            // std::cout << "Making 3D segment: "
            // 	      << " start(" << seg3d.start[0] << "," << seg3d.start[1] << "," << seg3d.start[2] << ") "
            // 	      << " end(" << seg3d.end[0] << "," << seg3d.end[1] << "," << seg3d.end[2] << ")"
            // 	      << std::endl;
            segments.emplace_back( std::move(seg3d) );
          }
        }
      }//end of j loop
    }//end of i loop
  }

}
