#ifndef __SELEXAMPLE_CXX__
#define __SELEXAMPLE_CXX__

#include "SelExample.h"

#include "InterTool_Util/InterImageUtils.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"

#include <sstream>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

namespace llcv {

  void SelExample::Configure (const larcv::PSet &pset) {
    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;

    LLCV_DEBUG() << "end" << std::endl;
  }

  void SelExample::Initialize() {
    _fout->cd();
    _outtree = new TTree("Example","");

    AttachRSEV(_outtree);
    
    return;
  }

  double SelExample::Select() {
    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "=======================" << std::endl;
    LLCV_DEBUG() << "(RSEV)=("<< Run() << "," << SubRun() << "," << Event() << "," << VertexID() << ")" << std::endl;    
    LLCV_DEBUG() << "VTX=(" 
		 << Data().Vertex()->X() << "," 
		 << Data().Vertex()->Y() << "," 
		 << Data().Vertex()->Z() << ")" << std::endl;

    //
    // get the image centered around the vertex
    //
    
    size_t cropx = 400;
    size_t cropy = 400;


    auto mat_v = Image().Image<cv::Mat>(kImageADC,cropx,cropy);
    auto img_v = Image().Image<larcv::Image2D>(kImageADC,cropx,cropy);
    
    std::vector<cv::Mat> mat3d_v;
    mat3d_v.reserve(3);

    for(const auto mat : mat_v) {
      auto img = larocv::Threshold(*mat,10,255);
      auto img3d = As8UC3(img);
      mat3d_v.emplace_back(std::move(img3d));
    }

    //
    // get the reconstructed tracks
    //
    
    auto track_v = Data().Tracks();

    for (const auto track : track_v) {

      for(size_t pid=0; pid< track->NumberTrajectoryPoints(); ++pid) {
	const auto& pt = track->LocationAtPoint(pid);
	for(size_t plane=0; plane<3; ++plane) {
	  int px_x = kINVALID_INT;
	  int px_y = kINVALID_INT;
	  ProjectImage2D(img_v.at(plane)->meta(),
			 pt.X(),pt.Y(),pt.Z(),
			 px_x, px_y);

	  if (px_x < 0) continue;
	  if (px_y < 0) continue;

	  if (px_x >= cropx) continue;
	  if (px_y >= cropy) continue;

	  LLCV_DEBUG() << "@plane=" << plane 
		    << " (" << px_x << "," << px_y << ")=" 
		    << img_v.at(plane)->pixel(px_y,px_x) << std::endl;
	  auto& mat3d = mat3d_v[plane];
	  mat3d.at<cv::Vec3b>(px_x,cropy-px_y-1) = {0,0,255};
	}
      }
    }

    
    for(size_t plane=0; plane<3; ++plane) {
      auto& img3d = mat3d_v[plane];
      std::stringstream ss; 
      ss << "plane_img_" << VertexID() << "_" << plane << ".png";
      cv::imwrite(ss.str(),img3d);
    }


    _outtree->Fill();
    
    LLCV_DEBUG() << "=======================" << std::endl;
    LLCV_DEBUG() << "end" << std::endl;
    return 0.0;
  }
  
  void SelExample::Finalize() {
    LLCV_DEBUG() << "start" << std::endl;

    _fout->cd();
    _outtree->Write();

    LLCV_DEBUG() << "end" << std::endl;
  }

}


#endif
