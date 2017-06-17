#include "Segment3DAlgo.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv
#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  Segment3DAlgo::Segment3DAlgo()
  : verbosity(0) {
  }

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
      if ( verbosity>0 )
	std::cout << __FILE__ << ":" << __LINE__ << " hits on p" << p
		  << " lo=" << hits_low.size()
		  << " hi=" << hits_high.size()
		  << " segs=" << segments2d.size()
		  << std::endl;
      if ( verbosity>1 ) {
	std::cout << "hit low list: " << std::endl;
	for (int i=0; i<(int)hits_low.size(); i++)
	  std::cout << " " << hits_low[i];
	std::cout << std::endl;
	std::cout << "hit high list: " << std::endl;
	for (int i=0; i<(int)hits_high.size(); i++)
	  std::cout << " " << hits_high[i];
	std::cout << std::endl;
      }
      plane_segments2d.emplace_back( std::move(segments2d) );
    }

    std::vector< Segment3D_t > segments;
    combine2Dinto3D( plane_segments2d, img_v, badch_v, hit_neighborhood, thresholds, 0.9, segments );

    return segments;
  }

  std::vector<int> Segment3DAlgo::findHits( const larcv::Image2D& img, const int row, const float& threshold, const int min_hit_width ) {
    std::vector<int> hits;

    // simple on-off hit finder with min with

    // if row out of bounds, return empty container
    if ( row<0 || row>(int)img.meta().rows() )
      return hits;
    
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
          seg2d.frac_w_charge = frac;
          seg2d.npix_w_charge = nrows_w_charge;
          seg2d.npix_tot = num_rows;
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
    int dcol  = abs(high_col-low_col);

    bool hascharge = false;
    nrows_w_charge = 0;
    num_rows = 0; // should read numsteps

    // we need to handle horizontal and vertical cases
    // we want to step using the slowest axis
    if ( drow>=dcol) {
      // vertical segments
      float slope = float(high_col-low_col)/float(high_row-low_row);
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
    else {
      // horizontal segments
      float slope = float(high_row-low_row)/float(high_col-low_col);
      for (int c=0; c<=dcol; c++) {
        int pixc = low_col + c;
        int pixr = low_row + slope*c;
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
      }//end of col loop
    }//end of if horizontal
  }

/*
  void Segment3DAlgo::checkSegmentChargeWithVariation( const std::vector<float>& fixedend, std::vector<float>& variable_end,
    const larcv::Image2D& img, const larcv::Image2D& badch, const float variation_dist, const std::vector<float>& thresholds, float& good_frac ) {
    // we vary the variable end in 3-space and calculate the good_frac.
    // we return the variable_end with the best good fraction and most in the current direction

    // make variation list
    std::vector< std::vector<float> > varend_list;
    varend_list.reserve(9);
    varend_list.push_back( variable_end ); // null variation

    // we basically are trying variations that amount to 1-pixel changes in the different project axes
    // tick direction: (1,0,0)
    // Y-plane: (0,0,1)
    // U-plane: (0, -sqrt(3)/2, 0.5)
    // V-plane: (0, +sqrt(3)/2, 0.5)
    const double sr3 = sqrt(3)/2.0;
    double varbases[8][3] = { { 1.0, 0.0, 0.0},
                              {-1.0, 0.0, 0.0},
                              { 0.0, 0.0, 1.0},
                              { 0.0, 0.0,-1.0},
                              { 0.0,-sr3, 0.5},
                              { 0.0, sr3,-0.5},
                              { 0.0,-sr3,-0.5},
                              { 0.0, sr3, 0.5} };
    // amounts to a 9 point star
    for (int b=0; b<8; b++) {
      std::vector<float> pos(3);
      for (int i=0; i<3; i++)
        pos[i] = variable_end[i] + variation_dist*varbases[b][i];
      varend_list.push_back( pos );
    }

    // get direction
    float segdir(3,0);
    float segnorm = 0.;
    for (int i=0; i<3; i++) {
      segdir[i] = varpos[i]-fixedend[i];
      segnorm += segdir[i]*segdir[i];
    }
    segnorm = sqrt(segnorm);
    for (int i=0; i<3; i++)
      segdir[i] /= segnorm;

    // now test
    // anchor image coordinates
    const larcv::ImageMeta& meta = img_v.front().meta();
    std::vector<int> fixedcoords = larcv::UBWireTool::getProjectedImagePixel( fixedend, meta, 3);

    std::vector< float > best_var_pos(3);
    float best_frac_good = 0;
    float best_cos = 0;
    for ( auto& varpos : varend_list ) {
      // get image coordinates
      std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( varpos, meta, 3 );
      for (size_t p=0; p<3; p++) {
        int lo_row, hi_row, lo_col, hi_col;
        if ( fixedcoords[0]<imgcoords[0] ) {
          lo_row  = fixedcoords[0];
          lo_col  = fixedcoords[p+1];
          hi_row  = imgcoords[0];
          hi_col  = imgcoords[p+1];
        }
        else {
          lo_row = imgcoords[0];
          lo_col = imgcoords[p+1];
          hi_row = fixedcoords[0];
          hi_col = fixedcoords[p+1];
        }
        int nrows, nrows_w_charge;
        checkSegmentCharge( img_v[p], badch_v[p], lo_row, lo_col, hi_row, lo_row, hit_neighborhood, thresholds[p], nrows_w_charge, nrows );
        float frac = 0.0;
        if ( nrows==0 ) {
          // not useable
          continue;
        }
        else {
          frac = float(nrows_w_charge)/float(nrows);
          // dir
          float varnorm = 0.;
          float vardir[3] = {0};
          for (int i=0; i<3; i++) {
            vardir[i] = varpos[i]-fixedend[i];
            varnorm = vardir[i]*vardir[i];
          }
          varnorm = sqrt(varnorm);
          float varcos = 0.;
          for (int i=0; i<3; i++) {
            varcos += vardir[i]*segdir[i]/varnorm;
          }

          if ( frac>best_frac_good ) {
            best_frac_good = frac;
            best_var_pos = varpos;
            best_cos = varcos;
          }
          else if ( frac==best_frac_good ) {
            if ( varcos>best_cos ) {
              varcos = best_cos;
              best_var_pos = varpos;
              best_frac_good = frac;
            }
          }
        }//else number of pixels along segment is non-zero
      }//end of plane
    }

  }
*/

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

        if ( verbosity>1 ) {
          std::cout << "      combo(" << i << "," << j << ") "
                    << "planes=(" << idx_a.plane << "," << idx_b.plane << ") "
                    << "frac_w_q=(" <<seg_a.frac_w_charge << "," << seg_b.frac_w_charge << ")"
                    << std::endl;
          std::cout << "        seg a=(" << seg_a.row_low << "," << seg_a.col_low << ") -> (" << seg_a.row_high << "," << seg_a.col_high << ")" << std::endl;
          std::cout << "        seg b=(" << seg_b.row_low << "," << seg_b.col_low << ") -> (" << seg_b.row_high << "," << seg_b.col_high << ")" << std::endl;
        }

        // check if out of the image
        if ( seg_a.row_low<0 || seg_b.row_low<0 || seg_a.row_high<0 || seg_b.row_high<0 )
          continue;

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
            std::cout << "        same high and low row " << std::endl;
        }
        else if ( same_high_row ) {
          // we have to adjust the bottom ends to match at the same row
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
            std::cout << "        same high (" << seg_a.row_high << ") rows. adjusted low= " << row_low
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
          float ftarget_wire_a = ( drow_a!=0 ) ? (dcol_a/drow_a)*(target_row-seg_a.row_low) + low_wire_a : high_wire_a;
          float dcol_b = (high_wire_b-low_wire_b);
          float drow_b = (seg_b.row_high-seg_b.row_low);
          float ftarget_wire_b = ( drow_b!=0 ) ? (dcol_b/drow_b)*(target_row-seg_b.row_low) + low_wire_b : high_wire_b;
          high_wire_a = ftarget_wire_a;
          high_wire_b = ftarget_wire_b;
          row_high = target_row;
          if ( verbosity>1 )  {
            std::cout << "        same low (" << seg_a.row_low << ") rows. adjusted high= " << row_high
              << " A(p=" << segidx[i].plane << "): (" << seg_a.row_high << "," << seg_a.col_high << ")->(" << row_high << "," << high_wire_a << ")"
              << " B(p=" << segidx[j].plane << "): (" << seg_b.row_high << "," << seg_b.col_high << ")->(" << row_high << "," << high_wire_b << ")"
              << std::endl;
          }
        }
        else {
          std::cout << "      should not get here." << std::endl;
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
          if ( verbosity>1 )
          std::cout << "        low hit does not cross." << std::endl;
          continue; // these wires don't intersect anywhere
        }

        // ok so both wires cross.
        if ( otherwire_high<0 || otherwire_high>=img_v.at(otherplane_high).meta().max_x() ) {
          if ( verbosity>1 )
          std::cout << "        high is out of bounds" << std::endl;
          continue;
        }
        if ( otherwire_low<0 || otherwire_low>=img_v.at(otherplane_low).meta().max_x() ) {
          if ( verbosity>1 )
            std::cout << "        low is out of bounds" << std::endl;
          continue;
        }
        // let's check the other plane for charge
        int othercol_high = img_v.at(otherplane_high).meta().col( otherwire_high );
        int othercol_low  = img_v.at(otherplane_low).meta().col( otherwire_low );

        // check the segment on the other plane
        int nrows_w_charge = 0;
        int nrows;
        checkSegmentCharge( img_v.at(otherplane_high), badch_v.at(otherplane_high), row_low, othercol_low, row_high, othercol_high,
          hit_width, thresholds.at(otherplane_high), nrows_w_charge, nrows );
        if ( verbosity>1 ) {
          // std::cout << "     2d combo(" << i << "," << j << ") "
          // 	  << " rowhi=" << row_high << " (p" << idx_a.plane << "," << high_wire_a << ") (p" << idx_b.plane << "," << high_wire_b << ") "
          // 	  << " rowlo=" << row_low  << " (p" << idx_a.plane << "," << low_wire_a  << ") (p" << idx_b.plane << "," << low_wire_b << ") "
          std::cout << "        otherplane=" << otherplane_high << " (" << row_low << "," << othercol_low << ") -> (" << row_high << "," << othercol_high << ")"
                    << " nrows=" << nrows << " frac=" << float(nrows_w_charge)/nrows
                    << std::endl;
        }

        // note that for plane V, the bipolar shape can integrate the charge to zero
        // so we will relax the good fraction for it, for tracks going along the V wire
        bool along_uv = false;
        if ( (otherplane_high==1 || otherplane_high==0 ) && abs(othercol_low-othercol_high)<=3 )
          along_uv = true;

        if ( (along_uv==true) || ( nrows>0 && along_uv==false && float(nrows_w_charge)/nrows > good_frac) ) {
          // make segment3d
          Segment3D_t seg3d;
          seg3d.start[1] = poszy_high[1];
          seg3d.start[2] = poszy_high[0];
          seg3d.start[0] = (img_v.at(otherplane_high).meta().pos_y( row_high )-3200.0)*(larutil::LArProperties::GetME()->DriftVelocity()*0.5);
          seg3d.end[1]   = poszy_low[1];
          seg3d.end[2]   = poszy_low[0];
          seg3d.end[0]   = (img_v.at(otherplane_low).meta().pos_y(  row_low )-3200.0)*(larutil::LArProperties::GetME()->DriftVelocity()*0.5);
          seg3d.plane_frac_w_charge.resize(3,0);
          seg3d.plane_frac_w_charge[ idx_a.plane ] = seg_a.frac_w_charge;
          seg3d.plane_frac_w_charge[ idx_b.plane ] = seg_b.frac_w_charge;
          seg3d.plane_frac_w_charge[ otherplane_high ] = float(nrows_w_charge)/float(nrows);

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
            if ( verbosity>1 ) {
              std::cout << "        Making 3D segment: "
                        << "q(" << seg3d.plane_frac_w_charge[0] << "," << seg3d.plane_frac_w_charge[1] << "," << seg3d.plane_frac_w_charge[2] << ")"
                	      << " start(" << seg3d.start[0] << "," << seg3d.start[1] << "," << seg3d.start[2] << ") "
                	      << " end(" << seg3d.end[0] << "," << seg3d.end[1] << "," << seg3d.end[2] << ")"
                 	      << std::endl;
            }
            segments.emplace_back( std::move(seg3d) );
          }
        }
      }//end of j loop
    }//end of i loop
  }

}
