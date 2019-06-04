#ifndef __SELTRACKDIR_CXX__
#define __SELTRACKDIR_CXX__

#include <numeric>

#include "SelTrackDir.h"
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TFitResult.h"

namespace llcv {

  void SelTrackDir::Configure (const larcv::PSet &pset) {
    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;

    fmax_hit_radius = pset.get<float>("MaxHitRadius");
    fplane          = pset.get<size_t>("Plane",2);
    fedge_length    = pset.get<float>("MaxEdgeLength",50);
    
    auto tradius = pset.get<float>("TruncatedRadius",3);
    ftsigma  = pset.get<float>("TruncatedSigma",0.5);

    ffit_exclusion = pset.get<float>("FitStartExclusion",8.0);

    _TruncMean.setRadius(tradius);

    LLCV_DEBUG() << "end" << std::endl;
  }

  void SelTrackDir::Initialize() {
    LLCV_DEBUG() << "start" << std::endl;
    
    fouttree = new TTree("trackdir","");
    AttachRSEV(fouttree);
    fouttree->Branch("trk_dedx_vv"              , &trk_dedx_vv);
    fouttree->Branch("trk_tdedx_vv"             , &trk_tdedx_vv);
    fouttree->Branch("trk_range_vv"             , &trk_range_vv);
    fouttree->Branch("trk_length_v"             , &trk_length_v);
    fouttree->Branch("trk_npts_v"               , &trk_npts_v);

    fouttree->Branch("trk_mean_dedx_v"          , &trk_mean_dedx_v);
    fouttree->Branch("trk_slope_dedx_v"         , &trk_slope_dedx_v);
    fouttree->Branch("trk_median_dedx_v"        , &trk_median_dedx_v);
    fouttree->Branch("trk_tmean_dedx_v"         , &trk_tmean_dedx_v);
    fouttree->Branch("trk_tslope_dedx_v"        , &trk_tslope_dedx_v);
    fouttree->Branch("trk_tmedian_dedx_v"       , &trk_tmedian_dedx_v);

    fouttree->Branch("trk_start_mean_dedx_v"    , &trk_start_mean_dedx_v);
    fouttree->Branch("trk_start_slope_dedx_v"   , &trk_start_slope_dedx_v);
    fouttree->Branch("trk_start_median_dedx_v"  , &trk_start_median_dedx_v);
    fouttree->Branch("trk_start_tmean_dedx_v"   , &trk_start_tmean_dedx_v);
    fouttree->Branch("trk_start_tslope_dedx_v"  , &trk_start_tslope_dedx_v);
    fouttree->Branch("trk_start_tmedian_dedx_v" , &trk_start_tmedian_dedx_v);

    fouttree->Branch("trk_middle_mean_dedx_v"   , &trk_middle_mean_dedx_v);
    fouttree->Branch("trk_middle_slope_dedx_v"  , &trk_middle_slope_dedx_v);
    fouttree->Branch("trk_middle_median_dedx_v" , &trk_middle_median_dedx_v);
    fouttree->Branch("trk_middle_tmean_dedx_v"  , &trk_middle_tmean_dedx_v);
    fouttree->Branch("trk_middle_tslope_dedx_v" , &trk_middle_tslope_dedx_v);
    fouttree->Branch("trk_middle_tmedian_dedx_v", &trk_middle_tmedian_dedx_v);
    
    fouttree->Branch("trk_end_mean_dedx_v"      , &trk_end_mean_dedx_v);
    fouttree->Branch("trk_end_slope_dedx_v"     , &trk_end_slope_dedx_v);
    fouttree->Branch("trk_end_median_dedx_v"    , &trk_end_median_dedx_v);
    fouttree->Branch("trk_end_tmean_dedx_v"     , &trk_end_tmean_dedx_v);
    fouttree->Branch("trk_end_tslope_dedx_v"    , &trk_end_tslope_dedx_v);
    fouttree->Branch("trk_end_tmedian_dedx_v"   , &trk_end_tmedian_dedx_v);
    
    fouttree->Branch("trk_forward_float_chi_v" , &trk_forward_float_chi_v);
    fouttree->Branch("trk_backward_float_chi_v", &trk_backward_float_chi_v);
    fouttree->Branch("trk_forward_float_A_v"   , &trk_forward_float_A_v);
    fouttree->Branch("trk_backward_float_A_v"  , &trk_backward_float_A_v);
    fouttree->Branch("trk_forward_float_d_v"   , &trk_forward_float_d_v);
    fouttree->Branch("trk_backward_float_d_v"  , &trk_backward_float_d_v);
    
    fouttree->Branch("trk_forward_fixed_p_chi_v" , &trk_forward_fixed_p_chi_v);
    fouttree->Branch("trk_backward_fixed_p_chi_v", &trk_backward_fixed_p_chi_v);

    fouttree->Branch("trk_forward_fixed_p_A_v"   , &trk_forward_fixed_p_A_v);
    fouttree->Branch("trk_backward_fixed_p_A_v"  , &trk_backward_fixed_p_A_v);

    fouttree->Branch("trk_forward_fixed_p_d_v"   , &trk_forward_fixed_p_d_v);
    fouttree->Branch("trk_backward_fixed_p_d_v"  , &trk_backward_fixed_p_d_v);
    
    fouttree->Branch("trk_forward_fixed_p_start_chi_v" , &trk_forward_fixed_p_start_chi_v);
    fouttree->Branch("trk_backward_fixed_p_start_chi_v", &trk_backward_fixed_p_start_chi_v);

    fouttree->Branch("trk_forward_fixed_p_middle_chi_v" , &trk_forward_fixed_p_middle_chi_v);
    fouttree->Branch("trk_backward_fixed_p_middle_chi_v", &trk_backward_fixed_p_middle_chi_v);

    fouttree->Branch("trk_forward_fixed_p_end_chi_v" , &trk_forward_fixed_p_end_chi_v);
    fouttree->Branch("trk_backward_fixed_p_end_chi_v", &trk_backward_fixed_p_end_chi_v);

    fouttree->Branch("trk_forward_fixed_m_chi_v" , &trk_forward_fixed_m_chi_v);
    fouttree->Branch("trk_backward_fixed_m_chi_v", &trk_backward_fixed_m_chi_v);

    fouttree->Branch("trk_forward_fixed_m_A_v"   , &trk_forward_fixed_m_A_v);
    fouttree->Branch("trk_backward_fixed_m_A_v"  , &trk_backward_fixed_m_A_v);

    fouttree->Branch("trk_forward_fixed_m_d_v"   , &trk_forward_fixed_m_d_v);
    fouttree->Branch("trk_backward_fixed_m_d_v"  , &trk_backward_fixed_m_d_v);

    fouttree->Branch("trk_forward_fixed_m_start_chi_v" , &trk_forward_fixed_m_start_chi_v);
    fouttree->Branch("trk_backward_fixed_m_start_chi_v", &trk_backward_fixed_m_start_chi_v);

    fouttree->Branch("trk_forward_fixed_m_middle_chi_v" , &trk_forward_fixed_m_middle_chi_v);
    fouttree->Branch("trk_backward_fixed_m_middle_chi_v", &trk_backward_fixed_m_middle_chi_v);

    fouttree->Branch("trk_forward_fixed_m_end_chi_v" , &trk_forward_fixed_m_end_chi_v);
    fouttree->Branch("trk_backward_fixed_m_end_chi_v", &trk_backward_fixed_m_end_chi_v);

    fouttree->Branch("trk_proton_hypo_chi_v" , &trk_proton_hypo_chi_v);
    fouttree->Branch("trk_proton_hypo_fit_chi_v" , &trk_proton_hypo_fit_chi_v);
    fouttree->Branch("trk_muon_hypo_chi_v"   , &trk_muon_hypo_chi_v);

    fouttree->Branch("trk_proton_hypo_start_chi_v" , &trk_proton_hypo_start_chi_v);
    fouttree->Branch("trk_proton_hypo_fit_start_chi_v" , &trk_proton_hypo_fit_start_chi_v);
    fouttree->Branch("trk_muon_hypo_start_chi_v"   , &trk_muon_hypo_start_chi_v);

    fouttree->Branch("trk_proton_hypo_middle_chi_v" , &trk_proton_hypo_middle_chi_v);
    fouttree->Branch("trk_proton_hypo_fit_middle_chi_v" , &trk_proton_hypo_fit_middle_chi_v);
    fouttree->Branch("trk_muon_hypo_middle_chi_v"   , &trk_muon_hypo_middle_chi_v);

    fouttree->Branch("trk_proton_hypo_end_chi_v" , &trk_proton_hypo_end_chi_v);
    fouttree->Branch("trk_proton_hypo_fit_end_chi_v" , &trk_proton_hypo_fit_end_chi_v);
    fouttree->Branch("trk_muon_hypo_end_chi_v"   , &trk_muon_hypo_end_chi_v);

    LLCV_DEBUG() << "end" << std::endl;
    return;
  }

  double SelTrackDir::Select() {
    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "=======================" << std::endl;

    ResetTree();

    const larlite::vertex& vtx = *(Data().Vertex());
    const auto& phit_v = Data().Hits();

    std::vector<larlite::hit> hit_v;
    hit_v.reserve(phit_v.size());
    for (auto const& phit : phit_v ) 
      hit_v.push_back(*phit);

    std::vector<int> hitmask_v(hit_v.size(),1);
    LLCV_DEBUG() << "number of hits: " << hit_v.size() << std::endl;
    LLCV_DEBUG() << "Number of tracks: "  << Data().Tracks().size() << std::endl;    

    const auto& trk_v = Data().Tracks();
    std::vector<std::vector<float>> dedx_track_per_plane;

    ResizeTree(trk_v.size());
    
    float stride = 0.5;

    for(size_t ithsort=0; ithsort<trk_v.size(); ++ithsort) { 
      _trackhitsort.clear();

      const auto& lltrack = *(trk_v[ithsort]);
      _trackhitsort.buildSortedHitList(vtx, lltrack, hit_v, fmax_hit_radius, hitmask_v);

      dedx_track_per_plane.clear();
      dedx_track_per_plane.resize(3);
      
      _trackhitsort.getPathBinneddEdx(stride,stride,dedx_track_per_plane);
      
      const auto& bincenter_xyz = _trackhitsort.getBinCentersXYZ(fplane); 
      const auto& dedx_track = dedx_track_per_plane.at(fplane);
      
      auto& trk_dedx_v = trk_dedx_vv[ithsort];
      trk_dedx_v.resize(dedx_track.size(),-1);

      static std::vector<float> dedx_v;
      static std::vector<float> range_v;
      static std::vector<float> rrange_v;
      dedx_v.clear();
      range_v.clear();
      rrange_v.clear();
      dedx_v.reserve(dedx_track.size());
      range_v.reserve(dedx_track.size());
      rrange_v.reserve(dedx_track.size());

      auto& length = trk_length_v[ithsort];
      auto& npts   = trk_npts_v[ithsort];

      npts = 0;

      // store the track
     
      auto& out_track = Data().MakeTrack();

      for (size_t ipt=0; ipt<dedx_track.size(); ipt++) {
	if (ipt >= bincenter_xyz.size()) continue;

	auto dedx = dedx_track.at(ipt);
	auto dx = stride*(float)(ipt+1);
	
	if (dedx==0) continue;
	
	trk_dedx_v[ipt] = dedx;

	dedx_v.push_back(dedx);
	range_v.push_back(dx);

	const auto& bincenter = bincenter_xyz.at(ipt);
	TVector3 xyz(bincenter.at(0),bincenter.at(1),bincenter.at(2));
	out_track.add_vertex(xyz);
	out_track.add_direction(xyz);
      }

      out_track.add_dqdx(std::vector<double>(dedx_track.begin(),dedx_track.end()));

      length = range_v.back();

      if (dedx_v.size()<3) continue;

      for(auto r : range_v)
	rrange_v.push_back(range_v.back() - r);

      npts = (int) range_v.size();

      static std::vector<float> tdedx_v;
      tdedx_v.clear();
      tdedx_v.reserve(range_v.size());

      _TruncMean.CalcTruncMeanProfile(range_v, dedx_v, tdedx_v, ftsigma);
      
      auto& trk_tdedx_v = trk_tdedx_vv[ithsort];
      auto& trk_range_v = trk_range_vv[ithsort];

      trk_tdedx_v = tdedx_v;
      trk_range_v = range_v;

      //
      // calculate overall 
      //
      
      auto& trk_mean_dedx = trk_mean_dedx_v[ithsort];
      auto& trk_slope_dedx = trk_slope_dedx_v[ithsort];
      auto& trk_median_dedx = trk_median_dedx_v[ithsort];

      auto& trk_tmean_dedx = trk_tmean_dedx_v[ithsort];
      auto& trk_tslope_dedx = trk_tslope_dedx_v[ithsort];
      auto& trk_tmedian_dedx = trk_tmedian_dedx_v[ithsort];
      
      trk_mean_dedx   = _TruncMean.Mean(dedx_v);
      trk_slope_dedx  = _TruncMean.Slope(range_v,dedx_v);
      trk_median_dedx = _TruncMean.Median(dedx_v);

      trk_tmean_dedx   = _TruncMean.Mean(tdedx_v);
      trk_tslope_dedx  = _TruncMean.Slope(range_v,tdedx_v);
      trk_tmedian_dedx = _TruncMean.Median(tdedx_v);
      
      //
      // start, middle, end 
      //
      
      size_t s0, s1;
      size_t m0, m1;
      size_t e0, e1;

      s0 = s1 = kINVALID_SIZE;
      m0 = m1 = kINVALID_SIZE;
      e0 = e1 = kINVALID_SIZE;

      // indexed range
      float range_sz      = (float)range_v.size();
      float range_sz_half = range_sz / 2.0;

      // indexed edge
      float fedge_sz      = fedge_length / stride;
      float fedge_sz_half = fedge_sz / 2.0;

      // physical range
      float range = range_sz * stride;

      LLCV_DEBUG() << "range_sz=" << range_sz << " range_sz_half=" << range_sz_half << std::endl;
      LLCV_DEBUG() << "fedge_sz=" << fedge_sz << " fedge_sz_half=" << fedge_sz_half << std::endl;
      LLCV_DEBUG() << "range   =" << range    << " fedge_length =" << fedge_length  << std::endl;

      bool case0 = false;

      if (range >= (3 * fedge_length) and range_sz >= (3 * fedge_sz))
	case0 = true;
      
      if (case0) {
	LLCV_DEBUG() << "case0" << std::endl;
	// 50 cm edges & middle
	s0 = 0;
	s1 = (size_t) (fedge_sz - 1);
	m0 = (size_t) (range_sz_half) - (size_t)(fedge_sz_half);
	m1 = (size_t) (range_sz_half) + (size_t)(fedge_sz_half) - 1;
	e1 = range_sz - 1;
	e0 = e1 - (size_t)(fedge_sz + 1);
      }
      else {
	LLCV_DEBUG() << "case1" << std::endl;
	// just chop it into 3 parts
	float range_sz_3 = range_sz / 3.0;
	s0 = 0;
	s1 = (size_t) range_sz_3 - 1;
	m0 = s1 + 1;
	m1 = m0 + (size_t) range_sz_3 - 1;
	e0 = m1 + 1;
	e1 = range_sz - 1;
      }

      LLCV_DEBUG() << " s0=" << s0 << " s1=" << s1 
		   << " m0=" << m0 << " m1=" << m1 
		   << " e0=" << e0 << " e1=" << e1 << std::endl;

      static std::vector<float> start_dedx_v;
      static std::vector<float> start_tdedx_v;
      static std::vector<float> start_range_v;
      static std::vector<float> start_rrange_v;
      start_dedx_v.clear();
      start_tdedx_v.clear();
      start_range_v.clear();
      start_rrange_v.clear();
      start_dedx_v.reserve(s1-s0);
      start_tdedx_v.reserve(s1-s0);
      start_range_v.reserve(s1-s0);
      start_rrange_v.reserve(s1-s0);

      static std::vector<float> middle_dedx_v;
      static std::vector<float> middle_tdedx_v;
      static std::vector<float> middle_range_v;
      static std::vector<float> middle_rrange_v;
      middle_dedx_v.clear();
      middle_tdedx_v.clear();
      middle_range_v.clear();
      middle_rrange_v.clear();
      middle_dedx_v.reserve(m1-m0);
      middle_tdedx_v.reserve(m1-m0);
      middle_range_v.reserve(m1-m0);
      middle_rrange_v.reserve(m1-m0);

      static std::vector<float> end_dedx_v;
      static std::vector<float> end_tdedx_v;
      static std::vector<float> end_range_v;
      static std::vector<float> end_rrange_v;
      end_dedx_v.clear();
      end_tdedx_v.clear();
      end_range_v.clear();
      end_rrange_v.clear();
      end_dedx_v.reserve(e1-e0);
      end_tdedx_v.reserve(e1-e0);
      end_range_v.reserve(e1-e0);
      end_rrange_v.reserve(e1-e0);

      for(size_t id=s0; id<=s1; ++id) {
	start_dedx_v.push_back(dedx_v[id]);
	start_tdedx_v.push_back(tdedx_v[id]);
	start_range_v.push_back(range_v[id]);
	start_rrange_v.push_back(range_v.back() - range_v[id]);
      }

      for(size_t id=m0; id<=m1; ++id) {
	middle_dedx_v.push_back(dedx_v[id]);
	middle_tdedx_v.push_back(tdedx_v[id]);
	middle_range_v.push_back(range_v[id]);
	middle_rrange_v.push_back(range_v.back() - range_v[id]);
      }
      
      for(size_t id=e0; id<=e1; ++id) {
	end_dedx_v.push_back(dedx_v[id]);
	end_tdedx_v.push_back(tdedx_v[id]);
	end_range_v.push_back(range_v[id]);
	end_rrange_v.push_back(range_v.back() - range_v[id]);
      }

      // start
      auto& trk_start_mean_dedx = trk_start_mean_dedx_v[ithsort];
      auto& trk_start_slope_dedx = trk_start_slope_dedx_v[ithsort];
      auto& trk_start_median_dedx = trk_start_median_dedx_v[ithsort];

      auto& trk_start_tmean_dedx = trk_start_tmean_dedx_v[ithsort];
      auto& trk_start_tslope_dedx = trk_start_tslope_dedx_v[ithsort];
      auto& trk_start_tmedian_dedx = trk_start_tmedian_dedx_v[ithsort];
      
      trk_start_mean_dedx   = _TruncMean.Mean(start_dedx_v);
      trk_start_slope_dedx  = _TruncMean.Slope(start_range_v,start_dedx_v);
      trk_start_median_dedx = _TruncMean.Median(start_dedx_v);

      trk_start_tmean_dedx   = _TruncMean.Mean(start_tdedx_v);
      trk_start_tslope_dedx  = _TruncMean.Slope(start_range_v,start_tdedx_v);
      trk_start_tmedian_dedx = _TruncMean.Median(start_tdedx_v);

      // middle
      auto& trk_middle_mean_dedx = trk_middle_mean_dedx_v[ithsort];
      auto& trk_middle_slope_dedx = trk_middle_slope_dedx_v[ithsort];
      auto& trk_middle_median_dedx = trk_middle_median_dedx_v[ithsort];

      auto& trk_middle_tmean_dedx = trk_middle_tmean_dedx_v[ithsort];
      auto& trk_middle_tslope_dedx = trk_middle_tslope_dedx_v[ithsort];
      auto& trk_middle_tmedian_dedx = trk_middle_tmedian_dedx_v[ithsort];
      
      trk_middle_mean_dedx   = _TruncMean.Mean(middle_dedx_v);
      trk_middle_slope_dedx  = _TruncMean.Slope(middle_range_v,middle_dedx_v);
      trk_middle_median_dedx = _TruncMean.Median(middle_dedx_v);

      trk_middle_tmean_dedx   = _TruncMean.Mean(middle_tdedx_v);
      trk_middle_tslope_dedx  = _TruncMean.Slope(middle_range_v,middle_tdedx_v);
      trk_middle_tmedian_dedx = _TruncMean.Median(middle_tdedx_v);

      // end
      auto& trk_end_mean_dedx = trk_end_mean_dedx_v[ithsort];
      auto& trk_end_slope_dedx = trk_end_slope_dedx_v[ithsort];
      auto& trk_end_median_dedx = trk_end_median_dedx_v[ithsort];

      auto& trk_end_tmean_dedx = trk_end_tmean_dedx_v[ithsort];
      auto& trk_end_tslope_dedx = trk_end_tslope_dedx_v[ithsort];
      auto& trk_end_tmedian_dedx = trk_end_tmedian_dedx_v[ithsort];
      
      trk_end_mean_dedx   = _TruncMean.Mean(end_dedx_v);
      trk_end_slope_dedx  = _TruncMean.Slope(end_range_v,end_dedx_v);
      trk_end_median_dedx = _TruncMean.Median(end_dedx_v);

      trk_end_tmean_dedx   = _TruncMean.Mean(end_tdedx_v);
      trk_end_tslope_dedx  = _TruncMean.Slope(end_range_v,end_tdedx_v);
      trk_end_tmedian_dedx = _TruncMean.Median(end_tdedx_v);


      std::vector<float> tdedx_rev_v;
      tdedx_rev_v.resize(tdedx_v.size());
      for(size_t rid=0; rid<tdedx_v.size(); ++rid)
	tdedx_rev_v[rid] = tdedx_v.at(tdedx_v.size() - rid - 1);
      
      //
      // Do the fit
      //
      TGraph tg_forward((Int_t)tdedx_v.size(),range_v.data(),tdedx_rev_v.data());
      TGraph tg_backward((Int_t)tdedx_v.size(),range_v.data(),tdedx_v.data());
      
      TF1 ftf_float("ftf_float",ftf_float_model,range_v.front(),range_v.back(),2);
      TF1 ftf_fixed_p("ftf_fixed_p",ftf_fixed_model_p,range_v.front(),range_v.back(),1);
      TF1 ftf_fixed_m("ftf_fixed_m",ftf_fixed_model_m,range_v.front(),range_v.back(),1);
      
      // Fit forward
      ftf_float.SetParameter(0,10);
      ftf_float.SetParameter(1,-0.42);
      ftf_fixed_p.SetParameter(0,10);
      ftf_fixed_m.SetParameter(0,10);
      
      TFitResultPtr ftf_forward_float_ptr = tg_forward.Fit(&ftf_float,"SQN0");
      TFitResultPtr ftf_forward_fixed_p_ptr = tg_forward.Fit(&ftf_fixed_p,"SQN0");
      TFitResultPtr ftf_forward_fixed_m_ptr = tg_forward.Fit(&ftf_fixed_m,"SQN0");

      auto& trk_forward_float_chi = trk_forward_float_chi_v[ithsort];
      auto& trk_forward_float_A = trk_forward_float_A_v[ithsort];
      auto& trk_forward_float_d = trk_forward_float_d_v[ithsort];

      auto& trk_forward_fixed_p_chi = trk_forward_fixed_p_chi_v[ithsort];
      auto& trk_forward_fixed_p_A = trk_forward_fixed_p_A_v[ithsort];
      auto& trk_forward_fixed_p_d = trk_forward_fixed_p_d_v[ithsort];

      auto& trk_forward_fixed_m_chi = trk_forward_fixed_m_chi_v[ithsort];
      auto& trk_forward_fixed_m_A = trk_forward_fixed_m_A_v[ithsort];
      auto& trk_forward_fixed_m_d = trk_forward_fixed_m_d_v[ithsort];
      
      trk_forward_float_chi = ftf_forward_float_ptr->Chi2();
      trk_forward_float_chi/= (float) ftf_forward_float_ptr->Ndf();
      trk_forward_float_A   = ftf_forward_float_ptr->Value(0);
      trk_forward_float_d   = ftf_forward_float_ptr->Value(1);

      trk_forward_fixed_p_chi = ftf_forward_fixed_p_ptr->Chi2();
      trk_forward_fixed_p_chi/= (float) ftf_forward_fixed_p_ptr->Ndf();
      trk_forward_fixed_p_A   = ftf_forward_fixed_p_ptr->Value(0);
      trk_forward_fixed_p_d   = -0.42;

      trk_forward_fixed_m_chi = ftf_forward_fixed_m_ptr->Chi2();
      trk_forward_fixed_m_chi/= (float) ftf_forward_fixed_m_ptr->Ndf();
      trk_forward_fixed_m_A   = ftf_forward_fixed_m_ptr->Value(0);
      trk_forward_fixed_m_d   = -0.37;
     
      // Fit backward
      ftf_float.SetParameter(0,10);
      ftf_float.SetParameter(1,-0.42);
      ftf_fixed_p.SetParameter(0,10);
      ftf_fixed_m.SetParameter(0,10);

      TFitResultPtr ftf_backward_float_ptr = tg_backward.Fit(&ftf_float,"SQN0");
      TFitResultPtr ftf_backward_fixed_p_ptr = tg_backward.Fit(&ftf_fixed_p,"SQN0");
      TFitResultPtr ftf_backward_fixed_m_ptr = tg_backward.Fit(&ftf_fixed_m,"SQN0");

      auto& trk_backward_float_chi = trk_backward_float_chi_v[ithsort];
      auto& trk_backward_float_A = trk_backward_float_A_v[ithsort];
      auto& trk_backward_float_d = trk_backward_float_d_v[ithsort];

      auto& trk_backward_fixed_p_chi = trk_backward_fixed_p_chi_v[ithsort];
      auto& trk_backward_fixed_p_A = trk_backward_fixed_p_A_v[ithsort];
      auto& trk_backward_fixed_p_d = trk_backward_fixed_p_d_v[ithsort];

      auto& trk_backward_fixed_m_chi = trk_backward_fixed_m_chi_v[ithsort];
      auto& trk_backward_fixed_m_A = trk_backward_fixed_m_A_v[ithsort];
      auto& trk_backward_fixed_m_d = trk_backward_fixed_m_d_v[ithsort];
      
      trk_backward_float_chi = ftf_backward_float_ptr->Chi2();
      trk_backward_float_chi/= (float) ftf_backward_float_ptr->Ndf();
      trk_backward_float_A   = ftf_backward_float_ptr->Value(0);
      trk_backward_float_d   = ftf_backward_float_ptr->Value(1);

      trk_backward_fixed_p_chi = ftf_backward_fixed_p_ptr->Chi2();
      trk_backward_fixed_p_chi/= (float) ftf_backward_fixed_p_ptr->Ndf();
      trk_backward_fixed_p_A   = ftf_backward_fixed_p_ptr->Value(0);
      trk_backward_fixed_p_d   = -0.42;

      trk_backward_fixed_m_chi = ftf_backward_fixed_m_ptr->Chi2();
      trk_backward_fixed_m_chi/= (float) ftf_backward_fixed_m_ptr->Ndf();
      trk_backward_fixed_m_A   = ftf_backward_fixed_m_ptr->Value(0);
      trk_backward_fixed_m_d   = -0.37;
      
      auto& trk_forward_fixed_p_start_chi   = trk_forward_fixed_p_start_chi_v[ithsort];
      auto& trk_backward_fixed_p_start_chi  = trk_backward_fixed_p_start_chi_v[ithsort];
      auto& trk_forward_fixed_p_middle_chi  = trk_forward_fixed_p_middle_chi_v[ithsort];
      auto& trk_backward_fixed_p_middle_chi = trk_backward_fixed_p_middle_chi_v[ithsort];
      auto& trk_forward_fixed_p_end_chi     = trk_forward_fixed_p_end_chi_v[ithsort];
      auto& trk_backward_fixed_p_end_chi     = trk_backward_fixed_p_end_chi_v[ithsort];

      trk_forward_fixed_p_start_chi   = RChi2WithFixedFit(range_v, tdedx_rev_v, s0, s1, trk_forward_fixed_p_A , trk_forward_fixed_p_d);
      trk_backward_fixed_p_start_chi  = RChi2WithFixedFit(range_v, tdedx_v    , s0, s1, trk_backward_fixed_p_A, trk_backward_fixed_p_d);
      trk_forward_fixed_p_middle_chi  = RChi2WithFixedFit(range_v, tdedx_rev_v, m0, m1, trk_forward_fixed_p_A , trk_forward_fixed_p_d);
      trk_backward_fixed_p_middle_chi = RChi2WithFixedFit(range_v, tdedx_v    , m0, m1, trk_backward_fixed_p_A, trk_backward_fixed_p_d);
      trk_forward_fixed_p_end_chi     = RChi2WithFixedFit(range_v, tdedx_rev_v, e0, e1, trk_forward_fixed_p_A , trk_forward_fixed_p_d);
      trk_backward_fixed_p_end_chi    = RChi2WithFixedFit(range_v, tdedx_v    , e0, e1, trk_backward_fixed_p_A, trk_backward_fixed_p_d);
      
      auto& trk_forward_fixed_m_start_chi   = trk_forward_fixed_m_start_chi_v[ithsort];
      auto& trk_backward_fixed_m_start_chi  = trk_backward_fixed_m_start_chi_v[ithsort];
      auto& trk_forward_fixed_m_middle_chi  = trk_forward_fixed_m_middle_chi_v[ithsort];
      auto& trk_backward_fixed_m_middle_chi = trk_backward_fixed_m_middle_chi_v[ithsort];
      auto& trk_forward_fixed_m_end_chi     = trk_forward_fixed_m_end_chi_v[ithsort];
      auto& trk_backward_fixed_m_end_chi    = trk_backward_fixed_m_end_chi_v[ithsort];
      
      trk_forward_fixed_m_start_chi   = RChi2WithFixedFit(range_v, tdedx_rev_v, s0, s1, trk_forward_fixed_m_A , trk_forward_fixed_m_d);
      trk_backward_fixed_m_start_chi  = RChi2WithFixedFit(range_v, tdedx_v    , s0, s1, trk_backward_fixed_m_A, trk_backward_fixed_m_d);
      trk_forward_fixed_m_middle_chi  = RChi2WithFixedFit(range_v, tdedx_rev_v, m0, m1, trk_forward_fixed_m_A , trk_forward_fixed_m_d);
      trk_backward_fixed_m_middle_chi = RChi2WithFixedFit(range_v, tdedx_v    , m0, m1, trk_backward_fixed_m_A, trk_backward_fixed_m_d);
      trk_forward_fixed_m_end_chi     = RChi2WithFixedFit(range_v, tdedx_rev_v, e0, e1, trk_forward_fixed_m_A , trk_forward_fixed_m_d);
      trk_backward_fixed_m_end_chi    = RChi2WithFixedFit(range_v, tdedx_v    , e0, e1, trk_backward_fixed_m_A, trk_backward_fixed_m_d);
      
      //
      // calculate the chi2 w.r.t theoretical
      //
      
      trk_proton_hypo_chi_v[ithsort]        = ProtonPIDChi2(tdedx_v,rrange_v);
      trk_proton_hypo_start_chi_v[ithsort]  = ProtonPIDChi2(start_tdedx_v,start_rrange_v);
      trk_proton_hypo_middle_chi_v[ithsort] = ProtonPIDChi2(middle_tdedx_v,middle_rrange_v);
      trk_proton_hypo_end_chi_v[ithsort]    = ProtonPIDChi2(end_tdedx_v,end_rrange_v);

      trk_proton_hypo_fit_chi_v[ithsort]        = ProtonFitPIDChi2(tdedx_v,rrange_v);
      trk_proton_hypo_fit_start_chi_v[ithsort]  = ProtonFitPIDChi2(start_tdedx_v,start_rrange_v);
      trk_proton_hypo_fit_middle_chi_v[ithsort] = ProtonFitPIDChi2(middle_tdedx_v,middle_rrange_v);
      trk_proton_hypo_fit_end_chi_v[ithsort]    = ProtonFitPIDChi2(end_tdedx_v,end_rrange_v);

      trk_muon_hypo_chi_v[ithsort]        = MuonPIDChi2(tdedx_v,rrange_v);
      trk_muon_hypo_start_chi_v[ithsort]  = MuonPIDChi2(start_tdedx_v,start_rrange_v);
      trk_muon_hypo_middle_chi_v[ithsort] = MuonPIDChi2(middle_tdedx_v,middle_rrange_v);
      trk_muon_hypo_end_chi_v[ithsort]    = MuonPIDChi2(end_tdedx_v,end_rrange_v);

      
    } // end this track
    
    fouttree->Fill();

    LLCV_DEBUG() << "=======================" << std::endl;
    LLCV_DEBUG() << "end" << std::endl;
    return 0.0;
  }
  
  void SelTrackDir::Finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    fouttree->Write();
    LLCV_DEBUG() << "end" << std::endl;
    return;
  }

  void SelTrackDir::ResizeTree(size_t ntracks) {
    
    trk_dedx_vv.resize(ntracks);
    trk_tdedx_vv.resize(ntracks);
    trk_range_vv.resize(ntracks);

    trk_length_v.resize(ntracks);
    trk_npts_v.resize(ntracks);

    trk_mean_dedx_v.resize(ntracks);
    trk_slope_dedx_v.resize(ntracks);
    trk_median_dedx_v.resize(ntracks);

    trk_tmean_dedx_v.resize(ntracks);
    trk_tslope_dedx_v.resize(ntracks);
    trk_tmedian_dedx_v.resize(ntracks);

    trk_start_mean_dedx_v.resize(ntracks);
    trk_start_slope_dedx_v.resize(ntracks);
    trk_start_median_dedx_v.resize(ntracks);

    trk_start_tmean_dedx_v.resize(ntracks);
    trk_start_tslope_dedx_v.resize(ntracks);
    trk_start_tmedian_dedx_v.resize(ntracks);

    trk_middle_mean_dedx_v.resize(ntracks);
    trk_middle_slope_dedx_v.resize(ntracks);
    trk_middle_median_dedx_v.resize(ntracks);

    trk_middle_tmean_dedx_v.resize(ntracks);
    trk_middle_tslope_dedx_v.resize(ntracks);
    trk_middle_tmedian_dedx_v.resize(ntracks);

    trk_end_mean_dedx_v.resize(ntracks);
    trk_end_slope_dedx_v.resize(ntracks);
    trk_end_median_dedx_v.resize(ntracks);

    trk_end_tmean_dedx_v.resize(ntracks);
    trk_end_tslope_dedx_v.resize(ntracks);
    trk_end_tmedian_dedx_v.resize(ntracks);

    trk_forward_float_chi_v.resize(ntracks);
    trk_backward_float_chi_v.resize(ntracks);
    
    trk_forward_float_A_v.resize(ntracks);
    trk_backward_float_A_v.resize(ntracks);
    
    trk_forward_float_d_v.resize(ntracks);
    trk_backward_float_d_v.resize(ntracks);

    trk_forward_fixed_p_chi_v.resize(ntracks);
    trk_backward_fixed_p_chi_v.resize(ntracks);
    
    trk_forward_fixed_p_A_v.resize(ntracks);
    trk_backward_fixed_p_A_v.resize(ntracks);
    
    trk_forward_fixed_p_d_v.resize(ntracks);
    trk_backward_fixed_p_d_v.resize(ntracks);

    trk_forward_fixed_p_start_chi_v.resize(ntracks);
    trk_backward_fixed_p_start_chi_v.resize(ntracks);
    trk_forward_fixed_p_middle_chi_v.resize(ntracks);
    trk_backward_fixed_p_middle_chi_v.resize(ntracks);
    trk_forward_fixed_p_end_chi_v.resize(ntracks);
    trk_backward_fixed_p_end_chi_v.resize(ntracks);

    trk_forward_fixed_m_chi_v.resize(ntracks);
    trk_backward_fixed_m_chi_v.resize(ntracks);
    
    trk_forward_fixed_m_A_v.resize(ntracks);
    trk_backward_fixed_m_A_v.resize(ntracks);
    
    trk_forward_fixed_m_d_v.resize(ntracks);
    trk_backward_fixed_m_d_v.resize(ntracks);

    trk_forward_fixed_m_start_chi_v.resize(ntracks);
    trk_backward_fixed_m_start_chi_v.resize(ntracks);
    trk_forward_fixed_m_middle_chi_v.resize(ntracks);
    trk_backward_fixed_m_middle_chi_v.resize(ntracks);
    trk_forward_fixed_m_end_chi_v.resize(ntracks);
    trk_backward_fixed_m_end_chi_v.resize(ntracks);

    trk_proton_hypo_chi_v.resize(ntracks);
    trk_proton_hypo_fit_chi_v.resize(ntracks);
    trk_muon_hypo_chi_v.resize(ntracks);

    trk_proton_hypo_start_chi_v.resize(ntracks);
    trk_proton_hypo_fit_start_chi_v.resize(ntracks);
    trk_muon_hypo_start_chi_v.resize(ntracks);

    trk_proton_hypo_middle_chi_v.resize(ntracks);
    trk_proton_hypo_fit_middle_chi_v.resize(ntracks);
    trk_muon_hypo_middle_chi_v.resize(ntracks);

    trk_proton_hypo_end_chi_v.resize(ntracks);
    trk_proton_hypo_fit_end_chi_v.resize(ntracks);
    trk_muon_hypo_end_chi_v.resize(ntracks);
  }

  void SelTrackDir::ResetTree() {

    trk_dedx_vv.clear();
    trk_tdedx_vv.clear();
    trk_range_vv.clear();

    trk_length_v.clear();
    trk_npts_v.clear();

    trk_mean_dedx_v.clear();
    trk_slope_dedx_v.clear();
    trk_median_dedx_v.clear();

    trk_tmean_dedx_v.clear();
    trk_tslope_dedx_v.clear();
    trk_tmedian_dedx_v.clear();

    trk_start_mean_dedx_v.clear();
    trk_start_slope_dedx_v.clear();
    trk_start_median_dedx_v.clear();

    trk_start_tmean_dedx_v.clear();
    trk_start_tslope_dedx_v.clear();
    trk_start_tmedian_dedx_v.clear();

    trk_middle_mean_dedx_v.clear();
    trk_middle_slope_dedx_v.clear();
    trk_middle_median_dedx_v.clear();

    trk_middle_tmean_dedx_v.clear();
    trk_middle_tslope_dedx_v.clear();
    trk_middle_tmedian_dedx_v.clear();

    trk_end_mean_dedx_v.clear();
    trk_end_slope_dedx_v.clear();
    trk_end_median_dedx_v.clear();

    trk_end_tmean_dedx_v.clear();
    trk_end_tslope_dedx_v.clear();
    trk_end_tmedian_dedx_v.clear();

    trk_forward_float_chi_v.clear();
    trk_backward_float_chi_v.clear();
    
    trk_forward_float_A_v.clear();
    trk_backward_float_A_v.clear();
    
    trk_forward_float_d_v.clear();
    trk_backward_float_d_v.clear();

    trk_forward_fixed_p_chi_v.clear();
    trk_backward_fixed_p_chi_v.clear();
    
    trk_forward_fixed_p_A_v.clear();
    trk_backward_fixed_p_A_v.clear();
    
    trk_forward_fixed_p_d_v.clear();
    trk_backward_fixed_p_d_v.clear();

    trk_forward_fixed_p_start_chi_v.clear();
    trk_backward_fixed_p_start_chi_v.clear();
    trk_forward_fixed_p_middle_chi_v.clear();
    trk_backward_fixed_p_middle_chi_v.clear();
    trk_forward_fixed_p_end_chi_v.clear();
    trk_backward_fixed_p_end_chi_v.clear();

    trk_forward_fixed_m_chi_v.clear();
    trk_backward_fixed_m_chi_v.clear();
    
    trk_forward_fixed_m_A_v.clear();
    trk_backward_fixed_m_A_v.clear();
    
    trk_forward_fixed_m_d_v.clear();
    trk_backward_fixed_m_d_v.clear();

    trk_forward_fixed_m_start_chi_v.clear();
    trk_backward_fixed_m_start_chi_v.clear();
    trk_forward_fixed_m_middle_chi_v.clear();
    trk_backward_fixed_m_middle_chi_v.clear();
    trk_forward_fixed_m_end_chi_v.clear();
    trk_backward_fixed_m_end_chi_v.clear();

    trk_proton_hypo_chi_v.clear();
    trk_proton_hypo_fit_chi_v.clear();
    trk_muon_hypo_chi_v.clear();

    trk_proton_hypo_start_chi_v.clear();
    trk_proton_hypo_fit_start_chi_v.clear();
    trk_muon_hypo_start_chi_v.clear();

    trk_proton_hypo_middle_chi_v.clear();
    trk_proton_hypo_fit_middle_chi_v.clear();
    trk_muon_hypo_middle_chi_v.clear();

    trk_proton_hypo_end_chi_v.clear();
    trk_proton_hypo_fit_end_chi_v.clear();
    trk_muon_hypo_end_chi_v.clear();
    
    return;
  }
  
  Double_t SelTrackDir::ftf_float_model(Double_t *x, Double_t *par) {
    static Double_t fitval = 0;
    fitval = par[0]*TMath::Power(x[0],par[1]);
    return fitval;
  }
  
  Double_t SelTrackDir::ftf_fixed_model_p(Double_t *x, Double_t *par) {
    static Double_t fitval = 0;
    fitval = par[0]*TMath::Power(x[0],-0.42);
    return fitval;
  }

  Double_t SelTrackDir::ftf_fixed_model_m(Double_t *x, Double_t *par) {
    static Double_t fitval = 0;
    fitval = par[0]*TMath::Power(x[0],-0.37);
    return fitval;
  }

  float SelTrackDir::RChi2WithFixedFit(const std::vector<float>& range_v,
				       const std::vector<float>& obs_v,
				       size_t start,
				       size_t end,
				       float A,
				       float d) {

    //rchi2 -- reduced chi2

    float chi2 = 0.0;

    Double_t range[1];
    Double_t par[2];
    par[0] = (Double_t) A;
    par[1] = (Double_t) d;

    assert(start <= end);
    assert(end   <  range_v.size());

    for(size_t id=start; id<end; ++id) {
      range[0] = (Double_t) range_v[id];
      float est = ftf_float_model(range,par);
      float local_chi2 = obs_v[id] - est;
      local_chi2 *= local_chi2;
      local_chi2 /= est;
      chi2 += local_chi2;
    }

    chi2 /= ((float) (end-start));
    
    return chi2;
  }


  float SelTrackDir::ProtonPIDChi2(const std::vector<float>& obs_v, const std::vector<float>& rr_v) {
    float chi2 = 0.0;
    
    assert (obs_v.size() == rr_v.size());
    
    for(size_t id=0; id<obs_v.size(); ++id) {
      float est = _dEdxCalculator.ProtonEstimate(rr_v[id]);
      float local_chi2 = obs_v[id] - est;
      local_chi2 *= local_chi2;
      local_chi2 /= est;
      chi2 += local_chi2;
    }
    
    chi2 /= (float)obs_v.size();

    return chi2;
  }

  float SelTrackDir::ProtonFitPIDChi2(const std::vector<float>& obs_v, const std::vector<float>& rr_v) {
    float chi2 = 0.0;
    
    assert (obs_v.size() == rr_v.size());
    
    for(size_t id=0; id<obs_v.size(); ++id) {
      float est = _dEdxCalculator.ProtonEstimateFit(rr_v[id]);
      float local_chi2 = obs_v[id] - est;
      local_chi2 *= local_chi2;
      local_chi2 /= est;
      chi2 += local_chi2;
    }
    
    chi2 /= (float)obs_v.size();

    return chi2;
  }

  
  float SelTrackDir::MuonPIDChi2(const std::vector<float>& obs_v, const std::vector<float>& rr_v) {
    float chi2 = 0.0;
    
    assert (obs_v.size() == rr_v.size());
    
    for(size_t id=0; id<obs_v.size(); ++id) {
      float est = _dEdxCalculator.MuonEstimate(rr_v[id]);
      float local_chi2 = obs_v[id] - est;
      local_chi2 *= local_chi2;
      local_chi2 /= est;
      chi2 += local_chi2;
    }
    
    chi2 /= (float)obs_v.size();
    
    return chi2;
  }

}


#endif
  
