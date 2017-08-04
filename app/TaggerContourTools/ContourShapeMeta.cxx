#include "ContourShapeMeta.h"

namespace larlitecv {

  ContourShapeMeta::ContourShapeMeta( const std::vector<cv::Point>& contour, const larcv::ImageMeta& meta )
    : std::vector<cv::Point>(contour),
      m_meta(meta),
      m_start( cv::Point(0,0) ),
      m_end( cv::Point(0,0) )  {
    _fill_linefit_members();
    _build_bbox();
  }

  void ContourShapeMeta::_fill_linefit_members() {
    cv::Vec4f out_array;

    cv::fitLine( *this, out_array, cv::DIST_L2, 0, 0.01, 0.01 );
    m_dir.resize(2,0);
    
    // loop through contour points and get min and max on projected line
    if ( out_array[1]!=0 ) {
      float norm  = sqrt(out_array[0]*out_array[0] + out_array[1]*out_array[1]);
      m_dir[0] = out_array[0]/norm;
      m_dir[1] = out_array[1]/norm;
    }
    else {
      // vertical
      float norm  = sqrt(out_array[0]*out_array[0] + out_array[1]*out_array[1]);
      m_dir[0] = out_array[0]/norm;
      m_dir[1] = 0;
    }
    
    cv::Point maxpt(0,0);
    cv::Point minpt(0,0);
    float mincos =  1e6;
    float maxcos = -1e6;
    for ( auto& pt : (*this) ) {
      float dx[2];
      dx[0] = pt.x-out_array[2];
      dx[1] = pt.y-out_array[3];
      float ptcos = 0.;
      for (int i=0; i<2; i++)
	ptcos += dx[i]*m_dir[i];
      if ( ptcos < mincos ) {
	minpt = pt;
	mincos = ptcos;
      }
      if ( ptcos > maxcos ) {
	maxpt = pt;
	maxcos = ptcos;
      }
    }
    m_start = minpt;
    m_end   = maxpt;

    // orient: start is at low y
    if ( m_start.y > m_end.y ) {
      cv::Point temp = m_start;
      m_start = m_end;
      m_end = temp;
      for (int i=0; i<2; i++)
	m_dir[i] *= -1.;
    }
    
  }
  
  void ContourShapeMeta::_build_bbox() {
    m_bbox = cv::boundingRect( *this );
  }
  
}
