#ifndef __SELMICHELSTUDY_CXX__
#define __SELMICHELSTUDY_CXX__

#include "SelMichelStudy.h"

#include "InterTool_Util/InterImageUtils.h"

#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"

#include <sstream>

namespace llcv {

  void SelMichelStudy::Configure (const larcv::PSet &pset) {
    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;

    LLCV_DEBUG() << "end" << std::endl;
  }

  double SelMichelStudy::Select() {
    LLCV_DEBUG() << "start" << std::endl;

    auto wire = Tree().Scalar<float>("wire");
    auto tick = Tree().Scalar<float>("tick");
    
    if (wire<0) return 0.0;
    if (tick<0) return 0.0;

    auto trk_mat_v = Image().Image<cv::Mat>(kImageTrack,-1,-1);
    auto shr_mat_v = Image().Image<cv::Mat>(kImageShower,-1,-1);

    const auto& pl2_trk_img = *(trk_mat_v[2]);
    const auto& pl2_shr_img = *(shr_mat_v[2]);

    tick /= 6;

    cv::Point_<float> pt(tick,wire);

    LLCV_INFO() << "michel @ " << pt << std::endl;

    int width = 200;
    int height = 200;

    auto small_trk_img = larocv::SmallImg(pl2_trk_img,pt,width,height);
    auto small_shr_img = larocv::SmallImg(pl2_shr_img,pt,width,height);

    small_trk_img = larocv::Threshold(small_trk_img,10,255);
    small_shr_img = larocv::Threshold(small_shr_img,10,255);

    small_trk_img = As8UC3(small_trk_img,2);
    small_shr_img = As8UC3(small_shr_img,0);

    auto small_img = small_trk_img + small_shr_img;

    std::stringstream ss;
    ss << "small_img_";
    ss << Run() << "_" << SubRun() << "_" << Event() << ".png";
    
    cv::imwrite(ss.str(),small_img);
    
    LLCV_DEBUG() << "end" << std::endl;
    return 0.0;
  }
  
  void SelMichelStudy::Finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    
    LLCV_DEBUG() << "end" << std::endl;
  }

}


#endif
