#include "ThruMuFoxTrot.h"

namespace larlitecv {

  // --------------------------------------------------------------------------
  // ThruMu Lead. Choose position most consistent with straight-line after
  // trying to remove space-charge effect.
  bool ThruMuFoxTrotLead::_chooseBestSegment_( const FoxStep& current, std::vector<Segment3D_t>& seg3d_v, const FoxTrack& past, const FoxTrotTrackerAlgoConfig& config, int& best_seg_idx ) {

    float closest_dist2realline = 1.0e9;
    best_seg_idx = -1;
    for ( int iseg=0; iseg<(int)seg3d_v.size(); iseg++ ) {
      Segment3D_t& seg = seg3d_v[iseg];

      // don't go backwards
      float segdist = 0.;
      float dirnorm = 0.;
      float segdir[3] = {0};
      for (int i=0; i<3; i++) {
	segdir[i] = seg.end[i]-seg.start[i];
	segdist += (seg.end[i]-seg.start[i])*(seg.end[i]-seg.start[i]);
	dirnorm += current.dir()[i]*current.dir()[i];
      }
      segdist = sqrt(segdist);
      dirnorm = sqrt(dirnorm);

      float cosdir = 0.;      
      for ( int i=0; i<3; i++ ){
	cosdir += current.dir()[i]*segdir[i]/(segdist*dirnorm);
      }


      // get original pos
      std::vector<double> segend(3,0);
      for (int i=0; i<3; i++)
	segend[i] = seg.end[i];
      segend[0] += m_shiftx;
      std::vector<double> doriginal_pos = m_sce.getOriginalPos( segend );

      double dist2realline = m_geoalgo.SqDist( real_line3d, doriginal_pos );
      
      std::cout << "seg3d #" << iseg << " dist2realline = " << dist2realline << " cosdir=" << cosdir << std::endl;      

      if ( cosdir>0.0 && dist2realline < closest_dist2realline ) {
	closest_dist2realline = dist2realline;
	best_seg_idx = iseg;
      }

    }//end of segment loop



    if ( best_seg_idx>=0 )
      return true;

    return false;
  }

  void ThruMuFoxTrotLead::setEndPoints( const std::vector<float>& start, const std::vector<float>& end, const float shiftx ) {

    m_shiftx = shiftx;
    
    std::vector<float> shift_start = start;
    shift_start[0] += m_shiftx;
    std::vector<float> shift_end   = end;
    shift_end[0]   += m_shiftx;
    
    image_line3d.Pt1( shift_start[0], shift_start[1], shift_start[2] );
    image_line3d.Pt2( shift_end[0],   shift_end[1],   shift_end[2] );

    std::vector<float> original_start = m_sce.getOriginalPos( shift_start );
    std::vector<float> original_end   = m_sce.getOriginalPos( shift_end );
    real_line3d.Pt1( original_start[0], original_start[1], original_start[2] );
    real_line3d.Pt2( original_end[0], original_end[1], original_end[2] );
  }

  void ThruMuFoxTrotLead::setEndPoints( const std::vector<double>& start, const std::vector<double>& end, const double shiftx ) {
    std::vector<float> fstart(3,0.0);
    std::vector<float> fend(3,0.0);
    for (int i=0; i<3; i++) {
      fstart[i] = start[i];
      fend[i] = end[i];
    }
    setEndPoints( fstart, fend, shiftx );
  }
  
  // --------------------------------------------------------------------------
  // THRUMU FOXTROT
  
  ThruMuFoxTrot::ThruMuFoxTrot( const ThruMuFoxTrotConfig& cfg )
    : m_config(cfg), m_tracker(cfg.foxtrotalgo_cfg) {
    m_tracker.setUserLead( &m_lead );
  }


  BMTrackCluster3D ThruMuFoxTrot::findThruMuTrack( const BoundarySpacePoint& pt_a, const BoundarySpacePoint& pt_b, 
						   const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v ) {
    BMTrackCluster3D bmtrack;
    return bmtrack;
  }
  
  FoxTrack ThruMuFoxTrot::findThruMuFoxTrack( const BoundarySpacePoint& pt_a, const BoundarySpacePoint& pt_b,
					      const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v ) {
    // we give the end points to the FoxTrot lead

    std::cout << "=========== THRUMU FOX TROT =======================" << std::endl;
    
    float shiftx = 0.;
    if ( pt_a.type()==larlitecv::kAnode ) {
      shiftx = -pt_a.pos()[0]; // should be zero
    }
    else if ( pt_b.type()==larlitecv::kAnode ) {
      shiftx = -pt_b.pos()[0]; // should be zero
    }
    else if ( pt_a.type()==larlitecv::kCathode ) {
      shiftx = 258.0 - pt_a.pos()[0];
    }
    else if ( pt_b.type()==larlitecv::kCathode ) {
      shiftx = 258.0 - pt_b.pos()[0];
    }
    std::cout << "Shift-x: " << shiftx << std::endl;
    
    m_lead.setEndPoints( pt_a.pos(), pt_b.pos(), shiftx );
    
    FoxTrack track = m_tracker.followTrack( img_v, badch_v, badch_v, pt_a );
    
    return track;
  }
  
}
