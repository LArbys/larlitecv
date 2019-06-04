#ifndef __SELTRACKDIR_H__
#define __SELTRACKDIR_H__

#include "InterTool_Core/InterSelBase.h"
#include "InterTool_Util/TruncMean.h"
#include "InterTool_Util/dEdxCalculator.h"
#include "TrackHitSorter/TrackHitSorter.h"


namespace llcv {
  
  class SelTrackDir : public InterSelBase { 

  public:

  SelTrackDir(std::string name="SelTrackDir") : 
    InterSelBase(name), 
    fouttree(nullptr)
    {}

    ~SelTrackDir(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize();
    double Select();
    void Finalize();
    void SetIsMC(int ismc)
    { bool ismc_b = (bool) ismc; _trackhitsort.SetIsMC(ismc_b); }
    
  private:

    TruncMean _TruncMean;
    dEdxCalculator _dEdxCalculator;
    larlitecv::TrackHitSorter _trackhitsort;

    TTree* fouttree;
    float fmax_hit_radius;
    size_t fplane;
    size_t fedge_length;
    float ftsigma;
    float ffit_exclusion;

    std::vector<std::vector<float> > trk_dedx_vv;
    std::vector<std::vector<float> > trk_tdedx_vv;
    std::vector<std::vector<float> > trk_range_vv;

    std::vector<float> trk_length_v;
    std::vector<int> trk_npts_v;

    //
    // overall
    //
    std::vector<float> trk_mean_dedx_v;
    std::vector<float> trk_slope_dedx_v;
    std::vector<float> trk_median_dedx_v;

    std::vector<float> trk_tmean_dedx_v;
    std::vector<float> trk_tslope_dedx_v;
    std::vector<float> trk_tmedian_dedx_v;

    //
    // start
    //
    std::vector<float> trk_start_mean_dedx_v;
    std::vector<float> trk_start_slope_dedx_v;
    std::vector<float> trk_start_median_dedx_v;

    std::vector<float> trk_start_tmean_dedx_v;
    std::vector<float> trk_start_tslope_dedx_v;
    std::vector<float> trk_start_tmedian_dedx_v;

    //
    // middle 
    //
    std::vector<float> trk_middle_mean_dedx_v;
    std::vector<float> trk_middle_slope_dedx_v;
    std::vector<float> trk_middle_median_dedx_v;

    std::vector<float> trk_middle_tmean_dedx_v;
    std::vector<float> trk_middle_tslope_dedx_v;
    std::vector<float> trk_middle_tmedian_dedx_v;
    
    //
    // end 
    //
    std::vector<float> trk_end_mean_dedx_v;
    std::vector<float> trk_end_slope_dedx_v;
    std::vector<float> trk_end_median_dedx_v;

    std::vector<float> trk_end_tmean_dedx_v;
    std::vector<float> trk_end_tslope_dedx_v;
    std::vector<float> trk_end_tmedian_dedx_v;


    //
    // dEdx = AR^(-d)
    //
    
    static Double_t ftf_float_model(Double_t *x, Double_t *par);
    static Double_t ftf_fixed_model_p(Double_t *x, Double_t *par);
    static Double_t ftf_fixed_model_m(Double_t *x, Double_t *par);

    std::vector<float> trk_forward_float_chi_v;
    std::vector<float> trk_backward_float_chi_v;

    std::vector<float> trk_forward_float_A_v;
    std::vector<float> trk_backward_float_A_v;

    std::vector<float> trk_forward_float_d_v;
    std::vector<float> trk_backward_float_d_v;

    std::vector<float> trk_forward_fixed_p_chi_v;
    std::vector<float> trk_backward_fixed_p_chi_v;

    std::vector<float> trk_forward_fixed_p_A_v;
    std::vector<float> trk_backward_fixed_p_A_v;

    std::vector<float> trk_forward_fixed_p_d_v;
    std::vector<float> trk_backward_fixed_p_d_v;

    std::vector<float> trk_forward_fixed_p_start_chi_v;
    std::vector<float> trk_backward_fixed_p_start_chi_v;

    std::vector<float> trk_forward_fixed_p_middle_chi_v;
    std::vector<float> trk_backward_fixed_p_middle_chi_v;

    std::vector<float> trk_forward_fixed_p_end_chi_v;
    std::vector<float> trk_backward_fixed_p_end_chi_v;    

    std::vector<float> trk_forward_fixed_m_chi_v;
    std::vector<float> trk_backward_fixed_m_chi_v;

    std::vector<float> trk_forward_fixed_m_A_v;
    std::vector<float> trk_backward_fixed_m_A_v;

    std::vector<float> trk_forward_fixed_m_d_v;
    std::vector<float> trk_backward_fixed_m_d_v;

    std::vector<float> trk_forward_fixed_m_start_chi_v;
    std::vector<float> trk_backward_fixed_m_start_chi_v;

    std::vector<float> trk_forward_fixed_m_middle_chi_v;
    std::vector<float> trk_backward_fixed_m_middle_chi_v;

    std::vector<float> trk_forward_fixed_m_end_chi_v;
    std::vector<float> trk_backward_fixed_m_end_chi_v;    

    
    //
    // chi2 against Theoretical
    //
    std::vector<float> trk_proton_hypo_chi_v;
    std::vector<float> trk_proton_hypo_fit_chi_v;
    std::vector<float> trk_muon_hypo_chi_v;

    std::vector<float> trk_proton_hypo_start_chi_v;
    std::vector<float> trk_proton_hypo_fit_start_chi_v;
    std::vector<float> trk_muon_hypo_start_chi_v;

    std::vector<float> trk_proton_hypo_middle_chi_v;
    std::vector<float> trk_proton_hypo_fit_middle_chi_v;
    std::vector<float> trk_muon_hypo_middle_chi_v;

    std::vector<float> trk_proton_hypo_end_chi_v;
    std::vector<float> trk_proton_hypo_fit_end_chi_v;
    std::vector<float> trk_muon_hypo_end_chi_v;
    
    void ResetTree();
    void ResizeTree(size_t ntracks);
    
    float RChi2WithFixedFit(const std::vector<float>& range_v,
			    const std::vector<float>& obs_v,
			    size_t start,
			    size_t end,
			    float A,
			    float d);

    float MuonPIDChi2(const std::vector<float>& obs_v, const std::vector<float>& rr_v);
    float ProtonPIDChi2(const std::vector<float>& obs_v, const std::vector<float>& rr_v);
    float ProtonFitPIDChi2(const std::vector<float>& obs_v, const std::vector<float>& rr_v);

  };

}

#endif

