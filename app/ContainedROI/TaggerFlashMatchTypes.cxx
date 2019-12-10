#include "TaggerFlashMatchTypes.h"

// LArCV
#include "DataFormat/DataFormatTypes.h"
#include "DataFormat/ROI.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// ROOT
#include "TVector.h"

namespace larlitecv {

  larcv::ROI TaggerFlashMatchData::MakeROI( const std::vector<larcv::Image2D>& img_v, const float bbox_pad_cm, const bool iscroi_candidate ) const {

    larcv::ROIType_t    roitype = larcv::kROIUnknown;
    larcv::ShapeType_t roishape = larcv::kShapeUnknown;

    if ( iscroi_candidate )
      roitype = larcv::kROIBNB;
    else
      roitype = larcv::kROICosmic;

    larcv::ROI roi( roitype, roishape );

    // we need the bounds in each wire plane
    float bb[3][2] = {0}; // extremes in each dimension
    float extrema[3][2][3] = { 0 }; // each extrema point stored here. (dim,min/max,xyz)
    for ( int v=0; v<3; v++) {
      bb[v][0] =  1.0e6;
      bb[v][1] = -1.0e6;
    }
  	
    for ( size_t i=0; i<m_track3d.NumberTrajectoryPoints(); i++ ) {
      const TVector3& xyz = m_track3d.LocationAtPoint(i);
      std::cout << "xyz: (" << xyz.X() << "," << xyz.Y() << "," << xyz.Z() << ")" << std::endl;
      for (int v=0; v<3; v++) {
	// minvalue
	if ( bb[v][0]>xyz[v] ) {
	  bb[v][0] = xyz[v];
	  for (int j=0; j<3; j++)
	    extrema[v][0][j] = xyz[j];
	}
	// maxvalue
	if ( bb[v][1]<xyz[v] ) {
	  bb[v][1] = xyz[v];
	  for (int j=0; j<3; j++)
	    extrema[v][1][j] = xyz[j];
	}
      }
    }

    // pad the box
    std::vector< std::vector<float> > corners(8);
    int icorner = 0;
    for (int i=0; i<2; i++) {
      for (int j=0; j<2; j++) {
	for (int k=0; k<2; k++) {
	  corners.at(icorner).resize(3,0);
	  corners[icorner][0] = bb[0][i] + (2*i-1)*bbox_pad_cm; // corner comes from different bboxes end
	  corners[icorner][1] = bb[1][j] + (2*j-1)*bbox_pad_cm;
	  corners[icorner][2] = bb[2][k] + (2*k-1)*bbox_pad_cm;
	  icorner++;
	}
      }
    }
    
    /// convert bounds into image coordinates
    const larcv::ImageMeta& meta = img_v.front().meta();
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    float min_tick = (bb[0][0]-bbox_pad_cm)/cm_per_tick + 3200.0;
    float max_tick = (bb[0][1]+bbox_pad_cm)/cm_per_tick + 3200.0;

    // pad
    min_tick -= 120.0;
    max_tick += 120.0;

    if ( ( min_tick < meta.min_y() || min_tick>meta.max_y() )
	 && ( max_tick < meta.min_y() || max_tick>meta.max_y() ) ) {
      std::cout << "[TaggerFlashMatchTypes::MakeROI] min (" << min_tick << ") and max (" << max_tick << ") tick out of range." << std::endl;
      std::cout << " verus meta: " << meta.dump() << std::endl;
      return roi;
    }

    if ( min_tick < meta.min_y() )
      min_tick = meta.min_y();
    if ( min_tick > meta.max_y() )
      min_tick = meta.max_y();

    if ( max_tick < meta.min_y() )
      max_tick = meta.min_y();
    if ( max_tick > meta.max_y() )
      max_tick = meta.max_y();

    int min_row = meta.row( max_tick );
    int max_row = meta.row( min_tick );

    std::vector<larcv::ImageMeta> bboxes;
    for (size_t p=0; p<3; p++) {

      // need to check each Y-Z.
      // find the wire-extrema for each
      float wire_extremes[2] = { 1.0e6, -1.0e6 };
      for ( auto const& corner : corners ) {
	Double_t xyz[3] = {0};
	for ( size_t i=0; i<3; i++ )
	  xyz[i] = corner[i];
	float wire = larutil::Geometry::GetME()->WireCoordinate(xyz,p);
	wire = (wire<0) ? 0 : wire;
	wire = (wire>=meta.max_x()) ? meta.max_x()-1.0 : wire;
	if ( wire<wire_extremes[0] )
	  wire_extremes[0] = wire;
	if ( wire>wire_extremes[1] )
	  wire_extremes[1] = wire;
      }

      //  enforce bounds
      if ( wire_extremes[0]<0 )
	wire_extremes[0] = 0;
      if ( p<2 && wire_extremes[1]>=2400 )
	wire_extremes[1] = 2399;
      if ( p==2 && wire_extremes[1]>=3456 )
	wire_extremes[1] = 3455;

      int col_min = meta.col( wire_extremes[0] );
      int col_max = meta.col( wire_extremes[1] );

      int ncols = abs(col_max-col_min);
      int nrows = abs(max_row-min_row);
      float width  = ncols*meta.pixel_width();
      float height = nrows*meta.pixel_height();

      larcv::ImageMeta planebb( width, height, nrows, ncols, wire_extremes[0], max_tick, (larcv::PlaneID_t)p );
      std::cout << planebb.dump();
      bboxes.emplace_back( std::move(planebb) );
    }

    roi.SetBB( bboxes );

    return roi;
  }
  

}
