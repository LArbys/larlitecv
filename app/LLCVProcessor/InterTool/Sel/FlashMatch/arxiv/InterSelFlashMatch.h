#ifndef __INTERSELFLASHMATCH_H__
#define __INTERSELFLASHMATCH_H__

#include "InterTool_Core/InterSelBase.h"
#include "FlashMatchInterface/GeneralFlashMatchAlgo.h"
#include "TrackHitSorter/TrackHitSorter.h"

// larlite
#include "DataFormat/vertex.h"
#include "DataFormat/track.h"

namespace llcv {
  
  class InterSelFlashMatch : public InterSelBase { 

  public:
    
  InterSelFlashMatch(std::string name="InterSelFlashMatch") : 
    InterSelBase(name)
      ,outtree(nullptr)
      ,genflashmatch(nullptr) 
    {}

    ~InterSelFlashMatch(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize();
    double Select();
    void Finalize();
    
    
  private:

    // -----------------------------
    // recorded trees
    TTree* outtree;


    // note: pmt values are opchannel-indexed
    
    float vtxpos_x;
    float vtxpos_y;
    float vtxpos_z;
    
    int number_tracks;
    int number_showers;

    // data flash
    int ndata_flashes;
    std::vector<float> data_totpe_v;
    std::vector<std::vector<float> > data_pe_vv;

    // track and track ID
    std::vector<std::vector<int> > proton_muon_pair_id_vv;
    std::vector<float> proton_muon_chi2_1mu1p_v;
    std::vector<float> proton_muon_chi2_shape_1mu1p_v;
    std::vector<float> proton_muon_hypo_totpe_v;
    std::vector<std::vector<float> > proton_muon_hypo_pe_vv;
    std::vector<int> proton_muon_best_data_flash_v;
    
    std::vector<std::vector<int> > muon_proton_pair_id_vv;
    std::vector<float> muon_proton_chi2_1mu1p_v;
    std::vector<float> muon_proton_chi2_shape_1mu1p_v;
    std::vector<float> muon_proton_hypo_totpe_v;
    std::vector<std::vector<float> > muon_proton_hypo_pe_vv;
    std::vector<int> muon_proton_best_data_flash_v;

    // track and shower ID
    std::vector<std::vector<int> > proton_shower_pair_id_vv;
    std::vector<float> proton_shower_chi2_1e1p_v;
    std::vector<float> proton_shower_chi2_shape_1e1p_v;
    std::vector<float> proton_shower_hypo_totpe_v;
    std::vector<std::vector<float> > proton_shower_hypo_pe_vv;
    std::vector<int> proton_shower_best_data_flash_v;

    void ResetVertex();

    // -----------------------------
    // running parameters
    larlitecv::GeneralFlashMatchAlgo* genflashmatch;
    float shower_correction_factor;
    float fmax_hit_radius;
    bool  isMC;
    bool  fSaveHistograms;

    flashana::QCluster_t build1e1pQCluster(const int protonid,
					   const larlite::vertex& vtx, 
					   const larlite::shower& shower, 
					   std::vector<larlitecv::TrackHitSorter>& dedxgen_v);

    flashana::QCluster_t build1mu1pQCluster(const int protonid,
					    const int muonid,
					    std::vector<larlitecv::TrackHitSorter>& dedxgen_v);

    void FillChi2(const std::vector<flashana::Flash_t>& dataflash_v,
		  const flashana::Flash_t& hypo,
		  float& best_chi2,
		  float& best_chi2_shape,
		  float& hypo_totpe,
		  std::vector<float>& hypo_pe,
		  int& data_flashidx);    

  };

}


#endif
