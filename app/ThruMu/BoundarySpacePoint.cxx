#include "BoundarySpacePoint.h"
#include <cmath>
#include <vector>
#include <sstream>
#include <string>

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv
#include "DataFormat/ImageMeta.h"
#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  float BoundarySpacePoint::dwall() const {
    float dwall = 0.0;
    float dy,dz;
    switch( type() ) {
    case larlitecv::kTop:
      dwall = 118.0-m_pos[1];
      break;
    case larlitecv::kBottom:
      dwall = 118.0+m_pos[1];
      break;
    case larlitecv::kUpstream:
      dwall = m_pos[2];
      break;
    case larlitecv::kDownstream:
      dwall = 1037-m_pos[2];
      break;
    case larlitecv::kAnode:
    case larlitecv::kCathode:
    case larlitecv::kImageEnd:
      dy = ( fabs(118.0-m_pos[1])<fabs(118.0+m_pos[1]) ) ? 118.0-m_pos[1] : 118.0+m_pos[1];
      dz = ( fabs(m_pos[2]) < fabs(1037-m_pos[2]) ) ? m_pos[2] : 1037-m_pos[2];
      dwall = ( fabs(dy)<fabs(dz ) ) ? dy : dz;
      break;
    default:
      std::runtime_error("BoundarySpacePoint::dwall[error] cannot calculate dwall for undefined boundary type");
      break;
    }
    return dwall;
  }

  void BoundarySpacePoint::setup( const larcv::ImageMeta& meta ) {
    if ( m_pos.size()==0 ) {
      m_pos.resize(3,0);
      m_dir.resize(3,0);

      // warning: unprotected use of meta
      
      float x = (meta.pos_y( front().row )-3200.0)*larutil::LArProperties::GetME()->DriftVelocity()*0.5;
      int crosses;
      std::vector<float> intersection;
      double triarea;
      std::vector< int > wids;
      for ( size_t p=0; p<size(); p++) {
	wids.push_back( meta.pos_x( at(p).col ) );
      }
      larcv::UBWireTool::wireIntersection( wids, intersection, triarea, crosses );
      m_pos[0] = x;
      m_pos[1] = intersection[1];
      m_pos[2] = intersection[2];
    }
  }

  int BoundarySpacePoint::tick( const larcv::ImageMeta& meta ) const {
    float tick = pos()[0]/(::larutil::LArProperties::GetME()->DriftVelocity()*0.5) + 3200.0;
    return (int)tick;
  }

  std::vector<int> BoundarySpacePoint::wires( const larcv::ImageMeta& meta ) const {
    std::vector<int> wires;
    Double_t xyz[3] = { pos()[0], pos()[1], pos()[2] };
    for (size_t p=0; p<3; p++) {
      wires.push_back( larutil::Geometry::GetME()->WireCoordinate( xyz, p ) );
    }
    return wires;
  }

  std::string BoundarySpacePoint::printImageCoords( const larcv::ImageMeta& meta ) const {
    int thetick = tick(meta);
    std::vector<int> thewires = wires(meta);
    std::stringstream ss;
    ss << "(" << thetick << ", " << thewires[0] << ", " << thewires[1] << ", " << thewires[2] << ")";
    return ss.str();
  }
}

// LEFT OVER CODE THAT TURNS PIXEL2D VECTORS INTO SPACEPOINTS. PROBABLY CRUFT.
/*
    // Start Point Information
    track3d.start_type = (larlitecv::BoundaryEnd_t) int(start_pt.front()->Intensity());
    track3d.row_start  = start_pt.front()->Y();
    track3d.tick_start = img_v.front().meta().pos_y( track3d.row_start );
    track3d.start_wire.resize(nplanes,0);
    track3d.start3D.resize(nplanes,0);
    for (int i=0; i<nplanes; i++) {
      track3d.start3D[i] = path.back().tyz[i];
      track3d.start_wire[i] = meta.pos_x( start_pt[i]->X() );
    }
    track3d.start3D[0] = (track3d.start3D[0]-3200)*cm_per_tick;

    // End Point Information
    if ( end_pt.size()!=img_v.size() ) {
      // No end specified
      track3d.end_type   = larlitecv::kUndefined;
      track3d.tick_end   = path.front().tyz.at(0);
      track3d.row_end    = meta.row( track3d.tick_end );
      track3d.end_wire.resize(nplanes,0);
      track3d.end3D.resize(nplanes,0);
      track3d.end3D[0] = (track3d.end3D[0]-3200)*cm_per_tick;
      for (int i=1; i<nplanes; i++)
        track3d.end3D[i] = path.front().tyz[i];
      Double_t xyz_end[3];
      for (int i=0; i<nplanes; i++)
        xyz_end[i] = track3d.end3D[i];
      for (int i=0; i<nplanes; i++) {
        float fwire = larutil::Geometry::GetME()->WireCoordinate(xyz_end,i);
        fwire = ( fwire<0 ) ? 0 : fwire;
        fwire = ( fwire>=meta.max_x() ) ? meta.max_x()-1.0 : fwire;
        track3d.end_wire[i] = (int)fwire;
      }
    }
    else {
      track3d.start_type = (larlitecv::BoundaryEnd_t) int(end_pt.front()->Intensity());
      track3d.row_end  = end_pt.front()->Y();
      track3d.tick_end = img_v.front().meta().pos_y( track3d.row_end );
      track3d.end_wire.resize(nplanes,0);
      track3d.end3D.resize(nplanes,0);
      for (int i=0; i<nplanes; i++) {
        track3d.end3D[i] = path.back().tyz[i];
        track3d.end_wire[i] = meta.pos_x( end_pt[i]->X() );
      }
      track3d.end3D[0] = (track3d.end3D[0]-3200)*cm_per_tick;
    }
    */
