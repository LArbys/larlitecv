#include "MCTrack2ImagePath.h"

// larlite
#include "LArUtil/LArProperties.h"

// larcv
#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  std::vector< std::vector<double> > mctrack2tyz( const larlite::mctrack& truthtrack,
						  const float trig_time,
						  bool returnxpos,
						  larlitecv::SpaceChargeMicroBooNE* psce ) {
    larlitecv::SpaceChargeMicroBooNE* sce = psce;
    if ( psce==NULL ) {
      // you want to avoid making one of these objects repeatedly
      sce = new larlitecv::SpaceChargeMicroBooNE;
    }

    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    // trig_time = 4050;
    
    // we loop through mcstep points in track
    // we project position in image, if inside image, we add to path

    std::vector< std::vector<double> > path_tyz_sce;
    
    for ( auto const& step : truthtrack ) {
      std::vector<float> pos(3);
      pos[0] = step.X();
      pos[1] = step.Y();
      pos[2] = step.Z();
      float t = step.T();

      // make space charge correction
      std::vector<double> pos_offset = sce->GetPosOffsets( pos[0], pos[1], pos[2] );
      std::vector<double> pos_sce(3);
      pos_sce[0] = pos[0]-pos_offset[0]+0.7;
      pos_sce[1] = pos[1]+pos_offset[1];
      pos_sce[2] = pos[2]+pos_offset[2];

      // set tick
      float tick = (t*1.0e-3 - (trig_time-4050.0))/0.5 + pos_sce[0]/cm_per_tick + 3200.0;

      if ( !returnxpos )
	pos_sce[0] = tick;
      else
	pos_sce[0] = (tick-3200.0)*cm_per_tick;
      path_tyz_sce.push_back( pos_sce );
    }

    return path_tyz_sce;
  }
  
  std::vector< std::vector<int> > tyzpath2imagepath( const std::vector< std::vector<double> >& path_tyz_sce, const std::vector<larcv::Image2D>& img_v ) {
    
    const larcv::ImageMeta& meta = img_v.front().meta();
    
    // now convert to image coordinates
    std::vector< std::vector<int> > imgpath;
    imgpath.reserve( path_tyz_sce.size() );
    for (auto const& pos : path_tyz_sce ) {
      std::vector<float> fpos(3);
      for (int i=0; i<3; i++)
	fpos[i] = pos[i];
      std::vector<int> crossing_imgcoords = larcv::UBWireTool::getProjectedImagePixel( fpos, meta, 3 );
      crossing_imgcoords[0] = (int)pos[0];
      imgpath.push_back( crossing_imgcoords );
    }
    
    return imgpath;
  }
  
  std::vector< std::vector<int> > mctrack2imagepath( const std::vector<larcv::Image2D>& img_v,
						     const larlite::mctrack& truthtrack,
						     const float trig_time,
						     larlitecv::SpaceChargeMicroBooNE* psce ) {
    std::vector< std::vector<double> > path_tyz_sce = mctrack2tyz( truthtrack, trig_time, false, psce );
    return tyzpath2imagepath( path_tyz_sce, img_v );
  }
}
