#include "RadialSegmentSearch.h"
#include <cstring>
#include <cmath>

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

#include "DataFormat/ImageMeta.h"
#include "UBWireTool/UBWireTool.h"

#include "RadialSegmentSearchTypes.h"

namespace larlitecv {

  std::vector< std::vector<larlitecv::RadialHit_t > > RadialSegmentSearch::findIntersectingChargeClusters( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
													   const std::vector<float>& pos3d, const float radius,
													   const std::vector<float>& thresholds ) {

    std::vector< std::vector<larcv::Pixel2D> > planerings = defineProjectedDiamondBoundary( pos3d, img_v, radius );
    std::vector< std::vector<larlitecv::RadialHit_t > > planehits;
    for (int p=0; p<(int)img_v.size(); p++) {
      std::vector< larlitecv::RadialHit_t > radhits = findIntersectingChargeClusters( img_v[p], badch_v[p], thresholds[p], planerings[p] );
      planehits.emplace_back( std::move(radhits) );
    }
    return planehits;
  }

  std::vector< larlitecv::RadialHit_t > RadialSegmentSearch::findIntersectingChargeClusters( const larcv::Image2D& img, const larcv::Image2D& badch,
											     const float threshold, const std::vector<larcv::Pixel2D>& pixelring ) {
    // Define a box using radius around the 3D position
    // project the wire range into the image
    // the projected box is also a box between the wire ranges and two ranges

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

    //std::cout << "starting ring index: " << idx_start << std::endl;
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
      if ( !inhit && idx+1!=idx_start) {
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
          loopedaround = false;
          //std::cout << "start hit at " << idx << " of " << pixelring.size() << " val=" << hit.maxval << " idx_start=" << idx_start << std::endl;
        }
        else {
          // do nothing
          //std::cout << " inhit: idx=" << idx << " val=" << pixelring[idx].Intensity() << std::endl;
        }
      }
      else if (inhit) {
        if ( pixelring[idx].Intensity()<threshold || idx+1==idx_start) {
          // end the hit
          if ( !loopedaround )
            hitlist.back().set_end( idx-1 );
          else
            hitlist.back().set_end( idx-1+(int)pixelring.size() );
          inhit = false;
          //std::cout << "end hit at " << idx << " of " << pixelring.size() << std::endl;
        }
        else {
          // update the max hit
          if ( !loopedaround )
            hitlist.back().update_max( idx, pixelring[idx].Intensity() );
          else
            hitlist.back().update_max( idx+(int)pixelring.size(), pixelring[idx].Intensity() );
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
    //std::cout << "return hitlist" << std::endl;
    return hitlist;
  }

  std::vector< std::vector< larcv::Pixel2D > > RadialSegmentSearch::defineProjectedDiamondBoundary( const std::vector<float>& center3d, const std::vector<larcv::Image2D>& img_v, const float radius ) {
    // we define a rhombus around the center3d point.
    // it is formed from the intersection of the U,V wires
    // we also define the box around the enclosing Y-wires

    const larcv::ImageMeta& meta = img_v.front().meta();
    Double_t upstreampt[3] = {center3d[0],center3d[1],center3d[2]};
    Double_t dnstreampt[3] = {center3d[0],center3d[1],center3d[2]};
    upstreampt[2] -= radius;
    dnstreampt[2] += radius;

    // we get the U and V wire bounds
    float wirebounds[3][2];
    for (size_t p=0; p<img_v.size(); p++) {
      wirebounds[p][0] = larutil::Geometry::GetME()->WireCoordinate( upstreampt, (larcv::PlaneID_t)p );
      wirebounds[p][1] = larutil::Geometry::GetME()->WireCoordinate( dnstreampt, (larcv::PlaneID_t)p );
      if ( wirebounds[p][0]<meta.min_x()  ) wirebounds[p][0] = 0;
      if ( wirebounds[p][0]>=meta.max_x() ) wirebounds[p][0] = meta.max_x()-1;
      if ( wirebounds[p][1]<meta.min_x()  ) wirebounds[p][1] = 0;            
      if ( wirebounds[p][1]>=meta.max_x() ) wirebounds[p][1] = meta.max_x()-1;

    }
    int colbounds[3][2] = {0};
    for (size_t p=0; p<3; p++) {
      colbounds[p][0] = meta.col( wirebounds[p][0] );
      colbounds[p][1] = meta.col( wirebounds[p][1] );
    }

    // get the time bounds
    float tickradius = radius/(larutil::LArProperties::GetME()->DriftVelocity()*0.5);
    float tick       = center3d[0]/(larutil::LArProperties::GetME()->DriftVelocity()*0.5)+3200.0;
    float tickbounds[2] = { tick-tickradius, tick+tickradius };
    int rowbounds[2] = {0};
    for (int i=0; i<2; i++) {
      if ( tickbounds[i]<=meta.min_y() ) tickbounds[i] = meta.min_y()+1;
      if ( tickbounds[i]>meta.max_y()  ) tickbounds[i] = meta.max_y();
      rowbounds[i] = meta.row( tickbounds[i] );
    }
    // rows and ticks are in reverse order, switch it
    int tmp = rowbounds[0];
    rowbounds[0] = rowbounds[1];
    rowbounds[1] = tmp;

    // with the wirebounds, we can define the pixel boundary in each image.
    // note that the boundary covers the spatial bounds for 2 sides of the diamond for u,v and the bottom of the diamond box
    // y will draw an actual box, but will overcover the diamond. assuming straight line segments, if it pierces y-boundary
    // we should be able to interpolate back to crossing

    std::vector< std::vector<larcv::Pixel2D> > pixelring;
    for ( size_t p=0; p<img_v.size(); p++ ) {
      std::vector<larcv::Pixel2D> planering;
      // make a box of pixels
      for (int c=colbounds[p][0]; c<=colbounds[p][1]; c++) {
        larcv::Pixel2D pix( c, rowbounds[1] ); // starts upper left, moves to upper right
        pix.Intensity( img_v[p].pixel(rowbounds[1],c) );
        planering.push_back( pix );
      }
      for (int r=rowbounds[1]; r>=rowbounds[0]; r--) {
        larcv::Pixel2D pix( colbounds[p][1], r );
        pix.Intensity( img_v[p].pixel(r,colbounds[p][1]) );
        planering.push_back( pix );
      }
      for (int c=colbounds[p][1]; c>=colbounds[p][0]; c--) {
        larcv::Pixel2D pix( c, rowbounds[0] );
        pix.Intensity( img_v[p].pixel(rowbounds[0],c) );
        planering.push_back( pix );
      }
      for (int r=rowbounds[0]; r<=rowbounds[1]; r++) {
        larcv::Pixel2D pix( colbounds[p][0], r );
        pix.Intensity( img_v[p].pixel(r,colbounds[p][0]) );
        planering.push_back( pix );
      }
      pixelring.emplace_back( std::move(planering) );
    }

    return pixelring;
  }

  std::vector< std::vector<larcv::Pixel2D> > RadialSegmentSearch::defineProjectedBoxBoundary( const std::vector<float>& center3d, const std::vector<larcv::Image2D>& img_v, const float radius ) {

    const larcv::ImageMeta& meta = img_v.front().meta();
    // first we build an order vector of pixels
    // they need to be unique, so we make an array to mark if we've used a pixel

    std::vector< std::vector<Double_t> > yzbox;
    for (int z=0; z<2; z++) {
      for (int y=0; y<2; y++) {
        std::vector<Double_t> yz(3,0.0);
        yz[1] = center3d[1] + (2*y-1)*radius;
        yz[2] = center3d[2] + (2*z-1)*radius;
        yzbox.push_back(yz);
      }
    }

    float colbounds[3][2] = {-1 };
    for (size_t p=0; p<3; p++) {
      for ( size_t i=0; i<yzbox.size(); i++ ) {
        float wire = larutil::Geometry::GetME()->WireCoordinate( yzbox[i], (int)p );
        wire = ( wire<0 ) ? 0 : wire;
        wire = (wire>=meta.max_x()) ? meta.max_x()-1 : wire;
        wire = (int)meta.col(wire);

        if ( colbounds[p][0]<0 || colbounds[p][0]>wire )
          colbounds[p][0] = wire;
        if ( colbounds[p][1]<0 || colbounds[p][1]<wire )
          colbounds[p][1] = wire;
      }
    }

    float tick = center3d[0]/(larutil::LArProperties::GetME()->DriftVelocity()*0.5)+3200.0;
    int row = meta.row( tick );
    int rowradius = radius/(larutil::LArProperties::GetME()->DriftVelocity()*0.5)/meta.pixel_height();

    std::vector< std::vector< larcv::Pixel2D > > planerings;
    for ( size_t p=0; p<3; p++ ) {

      // first define a bounding box
      std::vector<int> lowleft(2);
      lowleft[0] = colbounds[p][0];
      lowleft[1] = row-rowradius;
      lowleft[1] = ( lowleft[1]<0 ) ? 0 : lowleft[1];
      std::vector<int> upright(2);
      upright[0] = colbounds[p][1]+1;
      upright[1] = row+rowradius+1;
      upright[1] = ( upright[1]>=(int)meta.rows() ) ? (int)meta.rows()-1 : upright[1];

      std::vector< larcv::Pixel2D > pixelring;
      // make a box of pixels
      for (int c=colbounds[p][0]; c<=colbounds[p][1]; c++) {
        larcv::Pixel2D pix( c, upright[1] );
        pix.Intensity( img_v[p].pixel(upright[1],c) );
        pixelring.push_back( pix );
      }
      for (int r=upright[1]; r>=lowleft[1]; r--) {
        larcv::Pixel2D pix( colbounds[p][1], r );
        pix.Intensity( img_v[p].pixel(r,colbounds[p][1]) );
        pixelring.push_back( pix );
      }
      for (int c=colbounds[p][1]; c>=colbounds[p][0]; c--) {
        larcv::Pixel2D pix( c, lowleft[1] );
        pix.Intensity( img_v[p].pixel(lowleft[1],c) );
        pixelring.push_back( pix );
      }
      for (int r=lowleft[1]; r<=upright[1]; r++) {
        larcv::Pixel2D pix( colbounds[p][0], r );
        pix.Intensity( img_v[p].pixel(r,colbounds[p][0]) );
        pixelring.push_back( pix );
      }

      planerings.emplace_back( std::move(pixelring) );
    }

    return planerings;
  }

  std::vector< Segment2D_t > RadialSegmentSearch::make2Dsegments( const larcv::Image2D& img, const larcv::Image2D& badch, const std::vector<RadialHit_t>& hitlist,
								  const std::vector<float>& pos3d, const float threshold, const int min_hit_width,
								  const float frac_w_charges, const int verbosity ) {

    const int p = img.meta().plane();
    std::vector<int> img_coords = larcv::UBWireTool::getProjectedImagePixel( pos3d, img.meta(), 3 );
    Segment3DAlgo segalgo;

    std::vector< Segment2D_t > seg2d_v;
    for ( auto const& radhit : hitlist ) {
      Segment2D_t seg2d;
      //std::cout << "max idx=" << radhit.max_idx << std::endl;
      const larcv::Pixel2D& pix = radhit.pixlist.at( radhit.max_idx );
      if ( img_coords[0] < (int) pix.Y() ) {
        seg2d.row_high = (int)pix.Y();
        seg2d.col_high = (int)pix.X();
        seg2d.row_low  = img_coords[0];
        seg2d.col_low  = img_coords[p+1];
      }
      else {
        seg2d.row_low  = (int)pix.Y();
        seg2d.col_low  = (int)pix.X();
        seg2d.row_high = img_coords[0];
        seg2d.col_high = img_coords[p+1];
      }
      int rows_w_charge = 0;
      int num_rows = 0;
      segalgo.checkSegmentCharge( img, badch, seg2d.row_low, seg2d.col_low, seg2d.row_high, seg2d.col_high, min_hit_width, threshold, rows_w_charge, num_rows );
      float frac = (float)rows_w_charge/(float)num_rows;
      seg2d.frac_w_charge = frac;
      seg2d.npix_w_charge = rows_w_charge;
      seg2d.npix_tot = num_rows;
      // for debug
      if ( verbosity>2 ) {
        std::cout << "      radial segment: (" << seg2d.row_low << "," << seg2d.col_low << ") to (" << seg2d.row_high << "," << seg2d.col_high << ") "
                  << " frac charge=" << seg2d.frac_w_charge << " w/charge=" << seg2d.npix_w_charge << " nrows=" << seg2d.npix_tot << std::endl;
      }

      if ( frac > frac_w_charges ) {
        seg2d_v.emplace_back( std::move(seg2d) );
      }

    }

    return seg2d_v;
  }

  std::vector< Segment3D_t > RadialSegmentSearch::find3Dsegments( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
								  const std::vector<float>& pos3d, const float search_radius, const std::vector<float>& pixel_thresholds,
								  const int min_hit_width, const int hit_neighborhood, const float segment_frac_w_charge, int verbosity ) {
    // This function takes a 3D point, projects it into the plane, and returns possible 3D segments representing 3D line segments where there is charge

    // check arguments
    if ( img_v.size()!=badch_v.size() || img_v.size()!=pixel_thresholds.size() ) {
      throw std::runtime_error( "number of entries in img_v, badch_v, and/or pixel_thresholds does not match." );
    }
    if ( pos3d.size()!=3 )
      throw std::runtime_error( "expected 3D point for pos3d" );


    // first we find the hits around a 3d box of the 3D point
    float tick = pos3d[0]/(larutil::LArProperties::GetME()->DriftVelocity()*0.5) + 3200.0;
    int row = img_v.front().meta().row(tick);
    if ( verbosity>1 ) {
      std::cout << "  radialsegsearch -------------" << std::endl;
      std::vector<double> xyz(3);
      for (int i=0; i<3; i++)
        xyz[i] = pos3d[i];
      std::cout << "    current tick=" << tick << " row=" << row << " wires=(";
      for (int p=0; p<3; p++) {
        std::cout << larutil::Geometry::GetME()->WireCoordinate(xyz,p);
        if (p<2)
          std::cout << ", ";
      }
      std::cout << ")" << std::endl;
    }
    std::vector<int> origin_imgcoords = larcv::UBWireTool::getProjectedImagePixel( pos3d, img_v.front().meta(), 3 );

    std::vector< std::vector<Segment2D_t > > plane_seg2d_hi(img_v.size()); // segments that go up in row number (down in ticks)
    std::vector< std::vector<Segment2D_t > > plane_seg2d_lo(img_v.size()); // segments that go down in row number (up in ticks)

    m_planehits.clear();
    m_planehits = findIntersectingChargeClusters( img_v, badch_v, pos3d, search_radius, pixel_thresholds );

    for (size_t p=0; p<img_v.size(); p++) {

      // find potential charge hits along the surface of our box
      std::vector<RadialHit_t> radhits = m_planehits[p];
      if ( verbosity>1 ) {
        std::cout << "    radhits on p=" << p << ": " << radhits.size() << std::endl;
        if ( verbosity>2 ) {
          for (auto const& radhit : radhits ) {
            std::cout << "      radhit (r,c)=(" << radhit.pixlist[ radhit.max_idx ].Y() << "," << radhit.pixlist[ radhit.max_idx ].X() << ")"
        	      << " tick=" << img_v[p].meta().pos_y( radhit.pixlist[ radhit.max_idx ].Y() )
        	      << " val=" << radhit.maxval
        	      << std::endl;
          }
        }
      }

      // make 2D segments by simply drawing straight line from center point to the radial hits
      // we only keep those lines where a fraction of points contain charge
      std::vector<Segment2D_t> seg2d_v = make2Dsegments( img_v[p], badch_v[p], radhits, pos3d, pixel_thresholds[p],
							 hit_neighborhood, segment_frac_w_charge, verbosity );

      // we sort the segments into hi and lo
      for ( auto& seg2d : seg2d_v ) {
        if ( abs( seg2d.row_low-seg2d.row_high )>=2 ) {
          // moves across enough rows where we are confident in direction of segment
          if ( seg2d.row_low<row ) {
            plane_seg2d_lo[p].emplace_back( std::move(seg2d) );
          }
          else {
            plane_seg2d_hi[p].emplace_back( std::move(seg2d) );
          }
        }
        else {
          // horizontal segments need special treatment
          // boundary hit resolution probably not good to 1 pixel/row
          // we need to allow for many different combinations for 2D-3D matching
          if ( verbosity>1 )
            std::cout << "      horizontal segments" << std::endl;
          // segment goes into both hi and low
          if ( seg2d.row_low==row && abs(origin_imgcoords[p+1]-seg2d.col_low)<=1 ) {
            // the low point is the anchor of the segment
            // we create
            int drow =seg2d.row_high-seg2d.row_low;
            Segment2D_t segdn;
            segdn.row_high = seg2d.row_low;
            segdn.col_high = seg2d.col_low;
            segdn.row_low  = seg2d.row_low-drow;
            segdn.col_low  = seg2d.col_high;
            if ( verbosity> 1 ) {
              std::cout << "      low anchor: into hi vec: seg hi=(" << seg2d.row_high << "," << seg2d.col_high << ") lo=(" << seg2d.row_low << "," << seg2d.col_low << ")" << std::endl;
              std::cout << "      low anchor: into lo vec: seg hi=(" << segdn.row_high << "," << segdn.col_high << ") lo=(" << segdn.row_low << "," << segdn.col_low << ")" << std::endl;
            }
            plane_seg2d_hi[p].emplace_back( std::move(seg2d) );
            plane_seg2d_lo[p].emplace_back( std::move(segdn) );
          }
          else if ( seg2d.row_high==row && abs(origin_imgcoords[p+1]-seg2d.col_high)<1 ) {
            int drow =seg2d.row_high-seg2d.row_low;
            Segment2D_t segup;
            segup.row_low = seg2d.row_high;
            segup.col_low = seg2d.col_high;
            segup.row_high = seg2d.row_low+drow;
            segup.col_high = seg2d.col_low;
            if ( verbosity> 1 ) {
              std::cout << "      hi anchor: into hi vec: seg hi=(" << segup.row_high << "," << segup.col_high << ") lo=(" << segup.row_low << "," << segup.col_low << ")" << std::endl;
              std::cout << "      hi anchor: into lo vec: seg hi=(" << seg2d.row_high << "," << seg2d.col_high << ") lo=(" << seg2d.row_low << "," << seg2d.col_low << ")" << std::endl;
            }
            plane_seg2d_hi[p].emplace_back( std::move(segup) );
            plane_seg2d_lo[p].emplace_back( std::move(seg2d) );
          }
          else {
	    if ( verbosity> 1 ) {
	      std::cout << "      undefined case" << std::endl;
	    }
            // Segment2D_t segrightup  = seg2d;
            // segrightup.row_high++;
            // Segment2D_t segleftdn  = seg2d;
            // segleftdn.row_low--;
            // Segment2D_t segright = seg2d;
            // segleft.row_high++;
            // segdn.row_low--;
            // std::cout << "      sorting seg hi=(" << segup.row_high << "," << segup.col_high << ") lo=(" << segup.row_low << "," << segup.col_low << ")" << std::endl;
            // std::cout << "      sorting seg lo=(" << segdn.row_high << "," << segdn.col_high << ") lo=(" << segdn.row_low << "," << segdn.col_low << ")" << std::endl;
            //plane_seg2d_hi[p].emplace_back( std::move(segup) );
            //plane_seg2d_lo[p].emplace_back( std::move(segdn) );
            plane_seg2d_hi[p].push_back( std::move(seg2d) );
            plane_seg2d_lo[p].push_back( std::move(seg2d) );
          }
        }
      }
    }
    if ( verbosity>1 ) {
      std::cout << "    number of 2d segments: hi=("
                << plane_seg2d_hi[0].size()  << ","
                << plane_seg2d_hi[1].size() << ","
                << plane_seg2d_hi[2].size() << ") "
                << "lo=("
                << plane_seg2d_lo[0].size()  << ","
                << plane_seg2d_lo[1].size() << ","
                << plane_seg2d_lo[2].size() << ") "
                << std::endl;
    }

    // combine into 3D segments
    Segment3DAlgo algo;
    algo.setVerbosity( verbosity );
    std::vector< Segment3D_t > seg3d_v;
    algo.combine2Dinto3D( plane_seg2d_hi, img_v, badch_v, hit_neighborhood, pixel_thresholds, segment_frac_w_charge, seg3d_v );
    algo.combine2Dinto3D( plane_seg2d_lo, img_v, badch_v, hit_neighborhood, pixel_thresholds, segment_frac_w_charge, seg3d_v );

    return seg3d_v;
  }

}
