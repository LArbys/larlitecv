#include "Pixel2SpacePoint.h"

#include <sstream>

// larlite
#include "LArUtil/LArProperties.h"

// larcv
#include "UBWireTool/UBWireTool.h"

// larlitecv
#include "TaggerTypes/BoundaryEndPt.h"

namespace larlitecv {

  BoundarySpacePoint Pixel2SpacePoint( const std::vector<larcv::Pixel2D>& pixels, const BoundaryEnd_t endtype, const larcv::ImageMeta& meta ) {
    std::vector<BoundaryEndPt> endpt_v;
    std::vector<int> wire(pixels.size(),0);
    for (size_t p=0; p<pixels.size(); p++) {
      BoundaryEndPt pt( pixels[p].Y(), pixels[p].X(), endtype );
      wire[p] = meta.pos_x( pixels[p].X() );
      endpt_v.push_back(pt);
    }
    std::vector<float> poszy;
    int crosses = 0;
    double triarea;
    larcv::UBWireTool::wireIntersection( wire, poszy, triarea, crosses );
    if ( crosses==0 ) {
      std::stringstream msg;
      msg << __FILE__ << ":" << __LINE__ << ": Pixel combination leads to non-intersecting wires." << "\n";
      msg << "  Row=" << pixels.front().Y() << " Image Cols=(";
      for ( auto const& endpt : endpt_v ) {
	msg << endpt.col << " ";
      }
      msg << ")\n";
      msg << "  intersection (z,y): (" << poszy[0] << "," << poszy[1] << ") triarea=" << triarea << " crosses=" << crosses;
      msg << std::endl;
      
      throw std::runtime_error( msg.str() );
    }
    float x = ( meta.pos_y( pixels.front().Y() ) - 3200.0 )*(larutil::LArProperties::GetME()->DriftVelocity()*0.5);
    float y = poszy[1];
    float z = poszy[0];
    BoundarySpacePoint sp( endtype, std::move(endpt_v), x, y, z );
    return sp;
  }

}
