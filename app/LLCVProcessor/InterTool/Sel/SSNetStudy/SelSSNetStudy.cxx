#ifndef __SELSSNETSTUDY_CXX__
#define __SELSSNETSTUDY_CXX__

#include "SelSSNetStudy.h"

#include "InterTool_Util/InterImageUtils.h"

#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"

namespace llcv {

  void SelSSNetStudy::Configure (const larcv::PSet &pset) {
    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;
    ResetTree();
    LLCV_DEBUG() << "end" << std::endl;
  }

  void SelSSNetStudy::Initialize() {
    _fout->cd();
    _out_tree = new TTree("SelSSNetStudy","");

    AttachRSEV(_out_tree);
    _out_tree->Branch("ssnet_track_frac_vv", &_ssnet_track_frac_vv);
    _out_tree->Branch("ssnet_shower_frac_vv", &_ssnet_shower_frac_vv);
  }

  double SelSSNetStudy::Select() {
    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "=======================" << std::endl;
    llcv::logger::get_shared().set((msg::Level_t)2);

    //
    // just get the whole image
    //
    auto adc_mat_v  = Image().Image<cv::Mat>(kImageADC,-1,-1);
    auto trk_mat_v  = Image().Image<cv::Mat>(kImageTrack,-1,-1);
    auto shr_mat_v  = Image().Image<cv::Mat>(kImageShower,-1,-1);

    std::vector<cv::Mat> adc_mat_thresh_v;
    std::vector<cv::Mat> trk_mat_thresh_v;
    std::vector<cv::Mat> shr_mat_thresh_v;

    adc_mat_thresh_v.reserve(3);
    trk_mat_thresh_v.reserve(3);
    shr_mat_thresh_v.reserve(3);
    
    for(const auto mat : adc_mat_v) 
      adc_mat_thresh_v.emplace_back(larocv::Threshold(*mat,10,255));

    for(const auto mat : trk_mat_v) 
      trk_mat_thresh_v.emplace_back(larocv::Threshold(*mat,10,255));

    for(const auto mat : shr_mat_v) 
      shr_mat_thresh_v.emplace_back(larocv::Threshold(*mat,10,255));

    auto meta_v = Image().Image<larocv::ImageMeta>(kImageTrack,-1,-1);

    const auto& par_vv = Data().Particles();

    _ssnet_track_frac_vv.resize(par_vv.size());
    _ssnet_shower_frac_vv.resize(par_vv.size());

    for(size_t parid=0; parid < par_vv.size(); ++parid) {
      auto& ssnet_track_frac_v = _ssnet_track_frac_vv[parid];
      auto& ssnet_shower_frac_v = _ssnet_shower_frac_vv[parid];
      ssnet_track_frac_v.resize(3,-1);
      ssnet_shower_frac_v.resize(3,-1);
      
      const auto& par_v = par_vv[parid];
      
      for(size_t plane=0; plane<3; ++plane) {
	
	auto& ssnet_track_frac = ssnet_track_frac_v[plane];
	auto& ssnet_shower_frac = ssnet_shower_frac_v[plane];
	
	const auto& adc_mat = adc_mat_thresh_v[plane];
	const auto& trk_mat = trk_mat_thresh_v[plane];
	const auto& shr_mat = shr_mat_thresh_v[plane];
	const auto meta    = meta_v[plane];
	
	const auto& par = par_v[plane];
	if (par.empty()) continue;
	
	LLCV_CRITICAL() << "WARNING THIS MODULE IS NOT COMPATIBLE WITH 043018 UPDATE" << std::endl;
	throw llcv_err("die");

	auto ctor = AsContour(par,*meta);
	
	float trk_sz  = (float)larocv::CountNonZero(trk_mat,ctor,0);
	float shr_sz  = (float)larocv::CountNonZero(shr_mat,ctor,0);
	float ctor_sz = (float)larocv::CountNonZero(adc_mat,ctor,0);

	ssnet_track_frac  = trk_sz / ctor_sz;
	ssnet_shower_frac = shr_sz / ctor_sz;
      }
    }

    _out_tree->Fill();
    ResetTree();
    LLCV_DEBUG() << "=======================" << std::endl;
    LLCV_DEBUG() << "end" << std::endl;
    return 0.0;
  }
  

  void SelSSNetStudy::Finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    _fout->cd();
    _out_tree->Write();
    LLCV_DEBUG() << "end" << std::endl;
  }

  void SelSSNetStudy::ResetTree() {
    _ssnet_track_frac_vv.clear();
    _ssnet_shower_frac_vv.clear();
    return;
  }
}


#endif
