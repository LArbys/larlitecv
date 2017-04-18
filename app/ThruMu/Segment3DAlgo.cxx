#include "Segment3DAlgo.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv
#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  std::vector< Segment3D_t > Segment3DAlgo::find3DSegments( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
							  const int row_a, const int row_b, const std::vector<float>& thresholds, const int min_hit_width ) {

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
      std::vector< Segment2D_t > segments2d = make2DSegments( img_v.at(p), badch_v.at(p), lowrow, hits_low, highrow, hits_high, thresholds[p], min_hit_width, 0.9 );
      plane_segments2d.emplace_back( std::move(segments2d) );
    }

    std::vector< Segment3D_t > segments;
    combine2Dinto3D( plane_segments2d, img_v, badch_v, min_hit_width, thresholds, 0.9, segments );

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

  void Segment3DAlgo::checkSegmentCharge( const larcv::Image2D& img, const larcv::Image2D& badch, const int low_row, const int low_col, const int high_row, const int high_col, const int hit_width,
					  const float threshold, int& nrows_w_charge, int& num_rows ) {
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
      for (int dc=-hit_width; dc<=hit_width; dc++) {
	if ( pixc+dc<0 || pixc+dc>=(int)img.meta().cols() ) continue;
	if ( img.pixel(pixr,pixc+dc)>threshold || badch.pixel(pixr,pixc+dc)>0) {
	  hascharge = true;
	  break;
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

    for (int i=0; i<num_segs; i++) {
      for (int j=i+1; j<num_segs; j++) {
        // don't compare 2D segments on the same plane
        if ( segidx[i].plane==segidx[j].plane )
          continue;
        const SegID& idx_a = segidx[i];
        const SegID& idx_b = segidx[j];

        const Segment2D_t& seg_a = plane_segments2d.at(idx_a.plane).at(idx_a.idx);
        const Segment2D_t& seg_b = plane_segments2d.at(idx_b.plane).at(idx_b.idx);

        int high_wire_a = img_v.at(idx_a.plane).meta().pos_x( seg_a.col_high );
        int high_wire_b = img_v.at(idx_b.plane).meta().pos_x( seg_b.col_high );
        int low_wire_a  = img_v.at(idx_a.plane).meta().pos_x( seg_a.col_low );
        int low_wire_b  = img_v.at(idx_b.plane).meta().pos_x( seg_b.col_low );

        std::vector<float> poszy_high;
        int crosses_high = 0;
        int otherwire_high = 0;
        int otherplane_high = 0;
        larcv::UBWireTool::getMissingWireAndPlane( (int)idx_a.plane, high_wire_a, (int)idx_b.plane, high_wire_b, otherplane_high, otherwire_high, poszy_high, crosses_high );
        if ( crosses_high==0 )
          continue; // these wires don't intersect anywhere

        std::vector<float> poszy_low;
        int crosses_low = 0;
        int otherwire_low = 0;
        int otherplane_low = 0;
        larcv::UBWireTool::getMissingWireAndPlane( (int)idx_a.plane, low_wire_a, (int)idx_b.plane, low_wire_b, otherplane_low, otherwire_low, poszy_low, crosses_low );
        if ( crosses_low==0 )
          continue; // these wires don't intersect anywhere

        // ok so both wires cross.
        if ( otherwire_high<0 || otherwire_high>=img_v.at(otherplane_high).meta().max_x() ) continue;
        if ( otherwire_low<0 || otherwire_low>=img_v.at(otherplane_low).meta().max_x() ) continue;
        // let's check the other plane
        int othercol_high = img_v.at(otherplane_high).meta().col( otherwire_high );
        int othercol_low  = img_v.at(otherplane_low).meta().col( otherwire_low );

        int nrows_w_charge = 0;
        int nrows;
        checkSegmentCharge( img_v.at(otherplane_high), badch_v.at(otherplane_high), seg_a.row_low, othercol_low, seg_a.row_high, othercol_high, hit_width, thresholds.at(otherplane_high), nrows_w_charge, nrows );

        if ( nrows>0 && float(nrows_w_charge)/nrows > good_frac ) {
          // make segment3d

          Segment3D_t seg3d;
          seg3d.start[1] = poszy_high[1];
          seg3d.start[2] = poszy_high[0];
          seg3d.start[0] = (img_v.at(otherplane_high).meta().pos_y( seg_a.row_high )-3200.0)*(larutil::LArProperties::GetME()->DriftVelocity()*0.5);
          seg3d.end[1]   = poszy_low[1];
          seg3d.end[2]   = poszy_low[0];
          seg3d.end[0]   = (img_v.at(otherplane_low).meta().pos_y( seg_a.row_low )-3200.0)*(larutil::LArProperties::GetME()->DriftVelocity()*0.5);

	  // check it's not repeating
	  bool isduplicate = false;
	  for ( auto const& pastseg : segments ) {
	    float dist_hi = 0.;
	    float dist_lo = 0;
	    for (int i=0; i<3; i++) {
	      dist_hi = (seg3d.start[i]-pastseg.start[i])*(seg3d.start[i]-pastseg.start[i]);
	      dist_lo = (seg3d.end[i]-pastseg.end[i])*(seg3d.end[i]-pastseg.end[i]);	      
	    }
	    dist_hi = sqrt(dist_hi);
	    dist_lo = sqrt(dist_lo);
	    if ( dist_hi<3.0 && dist_lo<3.0 ) {
	      isduplicate = true;
	      break;
	    }
	  }

	  if ( !isduplicate ) {
	    std::cout << "Making 3D segment: "
		      << " start(" << seg3d.start[0] << "," << seg3d.start[1] << "," << seg3d.start[2] << ") "
		      << " end(" << seg3d.end[0] << "," << seg3d.end[1] << "," << seg3d.end[2] << ")"
		      << std::endl;
	    segments.emplace_back( std::move(seg3d) );
	  }
        }
      }//end of j loop
    }//end of i loop
  }

}
