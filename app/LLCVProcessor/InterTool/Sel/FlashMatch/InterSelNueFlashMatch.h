#ifndef __INTERSELNUEFLASHMATCH_H__
#define __INTERSELNUEFLASHMATCH_H__

#include "InterTool_Core/InterSelBase.h"
#include "FlashMatchInterface/GeneralFlashMatchAlgo.h"

// larlite
#include "DataFormat/track.h"

namespace llcv {
  
  class InterSelNueFlashMatch : public InterSelBase { 

  public:
    
  InterSelNueFlashMatch(std::string name="InterSelNueFlashMatch") : 
    InterSelBase(name)
      , outtree(nullptr)
      , genflashmatch(nullptr) 
    {}

    ~InterSelNueFlashMatch(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize();
    double Select();
    void Finalize();
    
  private:

    float fPXtoMEV;

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
    
    flashana::QCluster_t build1e1pQCluster(const larlite::track& trk1, 
					   const larlite::track& trk2);
    
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
