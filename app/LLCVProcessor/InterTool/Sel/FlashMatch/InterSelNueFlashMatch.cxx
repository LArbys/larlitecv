#ifndef __INTERSELNUEFLASHMATCH_CXX__
#define __INTERSELNUEFLASHMATCH_CXX__

#include "InterSelNueFlashMatch.h"

#include "TCanvas.h"
#include "TH1D.h"

#include <stdlib.h>

namespace llcv {

  void InterSelNueFlashMatch::Configure (const larcv::PSet &pset) {
    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;

    larcv::PSet genflash_pset = pset.get<larcv::PSet>("GeneralFlashMatchAlgo");
    larlitecv::GeneralFlashMatchAlgoConfig genflash_cfg = larlitecv::GeneralFlashMatchAlgoConfig::MakeConfigFromPSet( genflash_pset );
    genflashmatch = new larlitecv::GeneralFlashMatchAlgo( genflash_cfg );

    fPXtoMEV = pset.get<float>("PXtoMEV");
    
    LLCV_DEBUG() << "end" << std::endl;
  }

  void InterSelNueFlashMatch::Initialize() {
    if(!_fout) throw llcv_err("No output file?");

    _fout->cd();

    outtree = new TTree("ffmatch","");
    AttachRSEV(outtree);

    // rsev
    outtree->Branch("vtxpos_x",&vtxpos_x,"vtxpos_x/F");
    outtree->Branch("vtxpos_y",&vtxpos_y,"vtxpos_y/F");
    outtree->Branch("vtxpos_z",&vtxpos_z,"vtxpos_z/F");

    outtree->Branch("number_tracks",&number_tracks,"number_tracks/I");
    
    // data flash
    outtree->Branch("ndata_flashes", &ndata_flashes, "ndata_flashes/I");
    outtree->Branch("data_totpe_v", &data_totpe_v); // per flash
    outtree->Branch("data_pe_vv", &data_pe_vv); // per flash per channel
    
    // track and shower ID
    outtree->Branch("proton_shower_pair_vv",&proton_shower_pair_id_vv); // per pair it's a pair
    outtree->Branch("proton_shower_chi2_1e1p_v",&proton_shower_chi2_1e1p_v); // per pair
    outtree->Branch("proton_shower_chi2_shape_1e1p_v",&proton_shower_chi2_shape_1e1p_v);
    outtree->Branch("proton_shower_hypo_totpe_v",&proton_shower_hypo_totpe_v); // per pair
    outtree->Branch("proton_shower_hypo_pe_vv",&proton_shower_hypo_pe_vv); // per pair per channel
    outtree->Branch("proton_shower_best_data_flash_v",&proton_shower_best_data_flash_v); // per flash
  }

  double InterSelNueFlashMatch::Select() {
    LLCV_DEBUG() << "start" << std::endl;
    
    ResetVertex();
    
    const larlite::vertex& vtx = *(Data().Vertex());

    vtxpos_x = (float)vtx.X();
    vtxpos_y = (float)vtx.Y();
    vtxpos_z = (float)vtx.Z();
    
    LLCV_DEBUG() << "Number of tracks: "  << Data().Tracks().size() << std::endl;    
    
    // get tracks
    const auto& trk_v = Data().Tracks();
    
    number_tracks = (int)trk_v.size();
    
    //
    // track-shower hypothesis
    //
    static std::vector<int> pair_v(2,kINVALID_INT);
    
    std::vector<flashana::QCluster_t> qcluster_v;
    qcluster_v.reserve(trk_v.size());
    
    proton_shower_pair_id_vv.reserve(trk_v.size());

    for(size_t trkid1=0; trkid1<trk_v.size(); ++trkid1) {
      for(size_t trkid2=trkid1+1; trkid2<trk_v.size(); ++trkid2) {
	auto qcluster1  = build1e1pQCluster(*(trk_v[trkid1]), *(trk_v[trkid2]));
	auto qcluster2  = build1e1pQCluster(*(trk_v[trkid2]), *(trk_v[trkid1]));
	qcluster_v.emplace_back(std::move(qcluster1));
	qcluster_v.emplace_back(std::move(qcluster2));
	
	pair_v[0] = trkid1; 
	pair_v[1] = trkid2;
	proton_shower_pair_id_vv.push_back(pair_v);

	pair_v[0] = trkid2; 
	pair_v[1] = trkid1;
	proton_shower_pair_id_vv.push_back(pair_v);
      }
    }

    // Get data flashes
    const auto& opflash_ptr_v = Data().Flashes();

    std::vector<larlite::opflash> ev_opflash;
    ndata_flashes = (int)opflash_ptr_v.size();

    ev_opflash.reserve(ndata_flashes);
    for(const auto& opflash_ptr : opflash_ptr_v)
      ev_opflash.emplace_back(*opflash_ptr);

    data_totpe_v.resize(ndata_flashes,0.0);
    data_pe_vv.resize(ndata_flashes);
    for(auto& v : data_pe_vv) v.resize(32,kINVALID_FLOAT);

    std::vector<flashana::Flash_t> dataflash_v = genflashmatch->MakeDataFlashes(ev_opflash);

    for(size_t idata=0; idata<dataflash_v.size(); ++idata) {
      for (int ipmt=0; ipmt<32; ipmt++) {
	data_pe_vv[idata][ipmt] = dataflash_v[idata].pe_v[ipmt];
	data_totpe_v[idata] += dataflash_v[idata].pe_v[ipmt];
      }
    }
    
    //
    // make flash hypothesis for proton-shower pair
    //
    for(const auto& qcluster_1e1p : qcluster_v) {
      flashana::Flash_t hypo_1e1p = genflashmatch->GenerateUnfittedFlashHypothesis(qcluster_1e1p);

      proton_shower_chi2_1e1p_v.resize(proton_shower_chi2_1e1p_v.size()+1);
      auto& proton_shower_chi2_1e1p = proton_shower_chi2_1e1p_v.back();

      proton_shower_chi2_shape_1e1p_v.resize(proton_shower_chi2_shape_1e1p_v.size()+1);
      auto& proton_shower_chi2_shape_1e1p = proton_shower_chi2_shape_1e1p_v.back();

      proton_shower_hypo_totpe_v.resize(proton_shower_hypo_totpe_v.size()+1);
      auto& proton_shower_hypo_totpe = proton_shower_hypo_totpe_v.back();

      proton_shower_hypo_pe_vv.resize(proton_shower_hypo_pe_vv.size()+1);
      auto& proton_shower_hypo_pe_v = proton_shower_hypo_pe_vv.back();
      
      proton_shower_best_data_flash_v.resize(proton_shower_best_data_flash_v.size()+1);
      auto& proton_shower_best_data_flash = proton_shower_best_data_flash_v.back();

      proton_shower_hypo_pe_v.resize(32);

      FillChi2(dataflash_v,
	       hypo_1e1p,
	       proton_shower_chi2_1e1p,
	       proton_shower_chi2_shape_1e1p,
	       proton_shower_hypo_totpe,
	       proton_shower_hypo_pe_v,
	       proton_shower_best_data_flash);
    }

    outtree->Fill();
    
    LLCV_DEBUG() << "end" << std::endl;
    return 1;
  }

  void InterSelNueFlashMatch::FillChi2(const std::vector<flashana::Flash_t>& dataflash_v,
				       const flashana::Flash_t& hypo,
				       float& best_chi2,
				       float& best_chi2_shape,
				       float& hypo_totpe,
				       std::vector<float>& hypo_pe,
				       int& data_flashidx) {
    
    float maxpe_hypo = -1.0*kINVALID_FLOAT;    
    float petot_hypo = 0;

    hypo_totpe = 0;

    for (int ich=0; ich<32; ich++) {

      // fill array
      hypo_pe[ich] = hypo.pe_v[ich];
    
      if (hypo_pe[ich] < 1.0e-3) 
	hypo_pe[ich] = 1.0e-3;
    
      // total
      hypo_totpe += hypo_pe[ich];

      // max
      if (maxpe_hypo < hypo.pe_v[ich])

	maxpe_hypo = hypo.pe_v[ich];
      
    }
    maxpe_hypo /= petot_hypo;
    
    best_chi2 = -1*kINVALID_FLOAT; 

    // record which data flash best matched
    data_flashidx = -1;

    for (size_t idata=0; idata<dataflash_v.size(); idata++) {

      float chi2 = 0.;

      float totpe_data = 0.;

      for (int ipmt=0; ipmt<32; ipmt++)
	totpe_data += dataflash_v[idata].pe_v[ipmt];
      
      for (int ipmt=0; ipmt<32; ipmt++) {

	float pefrac_data = dataflash_v[idata].pe_v[ipmt] / totpe_data;
	float pefrac_hypo = hypo_pe[ipmt] / hypo_totpe;

	//float pe_data = pefrac_data * totpe_data;
	//float pefrac_data_err = std::sqrt(pe_data)/pe_data;

	float diff = pefrac_hypo - pefrac_data;
	
	// if ( pefrac_data_err>0 ) {
	//  chi2 += (diff*diff)/pefrac_data;
	// }
	// else if (pefrac_data_err==0.0 ){
	//  if ( pefrac_hypo>0 ) {
	//    chi2 += (diff*diff)/pefrac_hypo;
	//  }
	// }

	chi2 += (diff*diff)/pefrac_hypo;
	
      } // end pmt

      if ( (best_chi2<0) or (best_chi2>chi2) ) {
	best_chi2 = chi2;
	data_flashidx = idata;
      }
    }

    best_chi2_shape = best_chi2;
  }

  flashana::QCluster_t InterSelNueFlashMatch::build1e1pQCluster(const larlite::track& proton_track,
								const larlite::track& electron_track) {
    flashana::QCluster_t qinteraction;
    qinteraction.reserve(1000);
    
    std::vector<float> ly_proton_electron_v(2,-1);
    ly_proton_electron_v[0] = 19200;
    ly_proton_electron_v[1] = 20000;
    
    std::vector<const larlite::track*> proton_electron_v(2,nullptr);
    proton_electron_v[0] = &proton_track;
    proton_electron_v[1] = &electron_track;


    for(size_t tid=0; tid < proton_electron_v.size(); ++tid) {
      
      const auto& track = *(proton_electron_v[tid]);
      const auto ly     = ly_proton_electron_v[tid];
      
      bool has_yplane = false;
      size_t npts = track.NumberdQdx((larlite::geo::View_t)2);

      if ( npts > 0) {
	has_yplane = true;
      } 
      else {
	has_yplane = false;
	npts = track.NumberdQdx((larlite::geo::View_t)1);
      }

      LLCV_DEBUG() << "@tid=" << tid << " npts=" << npts << " yplane=" << (int)has_yplane << std::endl;

      const TVector3* loc_ptr = nullptr;
      float dedx = -1;

      for(size_t pid=0; pid < npts; ++pid) {
	
	if (has_yplane) {
	  loc_ptr = &(track.LocationAtPoint(pid)); 
	  dedx    = track.DQdxAtPoint(pid,(larlite::geo::View_t)2);
	}
	else {
	  loc_ptr = &(track.DirectionAtPoint(pid));
	  dedx    = track.DQdxAtPoint(pid,(larlite::geo::View_t)1);
	}

	if (dedx==0) continue;
      
	dedx /= fPXtoMEV;
      
	float numphotons = dedx * ly;

	flashana::QPoint_t qpt( loc_ptr->X(), loc_ptr->Y(), loc_ptr->Z(), numphotons );
	qinteraction.emplace_back( std::move(qpt) );

	LLCV_DEBUG() << "@tid=" << tid << ": (" << loc_ptr->X() << "," << loc_ptr->Y() << "," << loc_ptr->Z() << ") dedx=" << dedx << std::endl;
      }

    } // end track loop
    
    return qinteraction;
  }
  
  void InterSelNueFlashMatch::Finalize() {
    LLCV_DEBUG() << "start" << std::endl;

    _fout->cd();

    outtree->Write();

    LLCV_DEBUG() << "end" << std::endl;
  }

  void InterSelNueFlashMatch::ResetVertex() {
    
    vtxpos_x = kINVALID_FLOAT;
    vtxpos_y = kINVALID_FLOAT;
    vtxpos_z = kINVALID_FLOAT;

    number_tracks  = -1*kINVALID_INT;
    number_showers = -1*kINVALID_INT;
    
    // data flash
    ndata_flashes = -1.0*kINVALID_INT;
    data_totpe_v.clear();
    data_pe_vv.clear();

    // track and shower ID
    proton_shower_pair_id_vv.clear();
    proton_shower_chi2_1e1p_v.clear();
    proton_shower_chi2_shape_1e1p_v.clear();
    proton_shower_hypo_totpe_v.clear();
    proton_shower_hypo_pe_vv.clear();
    proton_shower_best_data_flash_v.clear();

  }

}


#endif
