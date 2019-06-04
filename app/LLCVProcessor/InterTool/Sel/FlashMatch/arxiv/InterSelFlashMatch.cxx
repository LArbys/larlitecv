#ifndef __INTERSELFLASHMATCH_CXX__
#define __INTERSELFLASHMATCH_CXX__

#include "InterSelFlashMatch.h"

#include "TCanvas.h"
#include "TH1D.h"

#include <stdlib.h>

namespace llcv {

  void InterSelFlashMatch::Configure (const larcv::PSet &pset) {
    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;
    larcv::PSet genflash_pset = pset.get<larcv::PSet>("GeneralFlashMatchAlgo");
    larlitecv::GeneralFlashMatchAlgoConfig genflash_cfg = larlitecv::GeneralFlashMatchAlgoConfig::MakeConfigFromPSet( genflash_pset );
    genflashmatch = new larlitecv::GeneralFlashMatchAlgo( genflash_cfg );
    shower_correction_factor = pset.get<float>("ShowerCorrectionFactor",1.0);
    isMC = pset.get<bool>("IsMC");
    fmax_hit_radius = pset.get<float>("MaxHitRadius");
    fSaveHistograms = pset.get<bool>("SaveFlashHistograms",false); // for debugging
    
    LLCV_DEBUG() << "end" << std::endl;
  }

  void InterSelFlashMatch::Initialize() {
    if(!_fout) throw llcv_err("No output file?");

    _fout->cd();

    outtree = new TTree("ffmatch","");
    AttachRSEV(outtree);

    // rsev
    outtree->Branch("vtxpos_x",&vtxpos_x,"vtxpos_x/F");
    outtree->Branch("vtxpos_y",&vtxpos_y,"vtxpos_y/F");
    outtree->Branch("vtxpos_z",&vtxpos_z,"vtxpos_z/F");

    outtree->Branch("number_tracks",&number_tracks,"number_tracks/I");
    outtree->Branch("number_showers",&number_showers,"number_showers/I");
    
    // data flash
    outtree->Branch("ndata_flashes", &ndata_flashes, "ndata_flashes/I");
    outtree->Branch("data_totpe_v", &data_totpe_v); // per flash
    outtree->Branch("data_pe_vv", &data_pe_vv); // per flash per channel

    // track and track ID
    outtree->Branch("proton_muon_pair_id_vv",&proton_muon_pair_id_vv);
    outtree->Branch("proton_muon_chi2_1mu1p_v",&proton_muon_chi2_1mu1p_v);
    outtree->Branch("proton_muon_chi2_shape_1mu1p_v",&proton_muon_chi2_shape_1mu1p_v);
    outtree->Branch("proton_muon_hypo_totpe_v",&proton_muon_hypo_totpe_v); // per pair
    outtree->Branch("proton_muon_hypo_pe_vv",&proton_muon_hypo_pe_vv); // per pair per channel
    outtree->Branch("proton_muon_best_data_flash_v",&proton_muon_best_data_flash_v); // per flash

    outtree->Branch("muon_proton_pair_id_vv",&muon_proton_pair_id_vv);
    outtree->Branch("muon_proton_chi2_1mu1p_v",&muon_proton_chi2_1mu1p_v);
    outtree->Branch("muon_proton_chi2_shape_1mu1p_v",&muon_proton_chi2_shape_1mu1p_v);
    outtree->Branch("muon_proton_hypo_totpe_v",&muon_proton_hypo_totpe_v); // per pair
    outtree->Branch("muon_proton_hypo_pe_vv",&muon_proton_hypo_pe_vv); // per pair per channel
    outtree->Branch("muon_proton_best_data_flash_v",&muon_proton_best_data_flash_v); // per flash

    // track and shower ID
    outtree->Branch("proton_shower_pair_vv",&proton_shower_pair_id_vv); // per pair it's a pair
    outtree->Branch("proton_shower_chi2_1e1p_v",&proton_shower_chi2_1e1p_v); // per pair
    outtree->Branch("proton_shower_chi2_shape_1e1p_v",&proton_shower_chi2_shape_1e1p_v);
    outtree->Branch("proton_shower_hypo_totpe_v",&proton_shower_hypo_totpe_v); // per pair
    outtree->Branch("proton_shower_hypo_pe_vv",&proton_shower_hypo_pe_vv); // per pair per channel
    outtree->Branch("proton_shower_best_data_flash_v",&proton_shower_best_data_flash_v); // per flash
    
  }

  double InterSelFlashMatch::Select() {
    LLCV_DEBUG() << "start" << std::endl;
    
    ResetVertex();
    
    const larlite::vertex& vtx = *(Data().Vertex());

    vtxpos_x = (float)vtx.X();
    vtxpos_y = (float)vtx.Y();
    vtxpos_z = (float)vtx.Z();

    // GET DATA
    // --------
    // get hits
    const auto& phit_v = Data().Hits();

    std::vector< larlite::hit > hit_v;
    hit_v.reserve(phit_v.size());
    for (auto const& phit : phit_v ) hit_v.push_back(*phit);

    std::vector<int> hitmask_v(hit_v.size(),1);
    LLCV_DEBUG() << "number of hits: " << hit_v.size() << std::endl;
    
    LLCV_DEBUG() << "Number of showers: " << Data().Showers().size() << std::endl;
    LLCV_DEBUG() << "Number of tracks: "  << Data().Tracks().size() << std::endl;    

    // get shower
    const auto& shr_v = Data().Showers();

    // get tracks
    const auto& trk_v = Data().Tracks();

    number_tracks = (int)trk_v.size();
    number_showers= (int)shr_v.size();

    // this class will generated de/dx per 3d position along the track. will use to set MeV deposited at 3d pos.
    std::vector<larlitecv::TrackHitSorter> dedxgen_v(trk_v.size());

    for(int ithsort=0; ithsort<(int)trk_v.size(); ++ithsort) { 
      const auto& lltrack = *(trk_v[ithsort]);
      dedxgen_v[ithsort].buildSortedHitList(vtx, lltrack, hit_v, fmax_hit_radius, hitmask_v);
    }

    // ok, we now have dedx 3d tracks for all tracks and shower
    // we want to build 1mu1p and 1e1p hypothesis

    //
    // track-track hypothesis
    //
    std::vector<flashana::QCluster_t> qcluster_proton_muon_v;
    qcluster_proton_muon_v.reserve(trk_v.size()*2);
    
    std::vector<flashana::QCluster_t> qcluster_muon_proton_v;
    qcluster_muon_proton_v.reserve(trk_v.size()*2);

    proton_muon_pair_id_vv.reserve(trk_v.size()*2);
    muon_proton_pair_id_vv.reserve(trk_v.size()*2);

    static std::vector<int> pair_v(2,kINVALID_INT);

    for(size_t trkid1=0; trkid1<trk_v.size(); ++trkid1) {
      for(size_t trkid2=trkid1+1; trkid2<trk_v.size(); ++trkid2) {
	flashana::QCluster_t qcluster_proton_muon = build1mu1pQCluster(trkid1, trkid2, dedxgen_v);
	flashana::QCluster_t qcluster_muon_proton = build1mu1pQCluster(trkid2, trkid1, dedxgen_v);

	qcluster_proton_muon_v.emplace_back(std::move(qcluster_proton_muon));
	qcluster_muon_proton_v.emplace_back(std::move(qcluster_muon_proton));

	pair_v[0] = trkid1;
	pair_v[1] = trkid2;
	proton_muon_pair_id_vv.push_back(pair_v);

	pair_v[0] = trkid2;
	pair_v[1] = trkid1;
	muon_proton_pair_id_vv.push_back(pair_v);
      }
    }

    std::vector<flashana::QCluster_t> qcluster_proton_shower_v;
    qcluster_proton_shower_v.reserve(trk_v.size() + shr_v.size());

    proton_shower_pair_id_vv.reserve(trk_v.size() + shr_v.size());

    for(size_t trkid1=0; trkid1<trk_v.size(); ++trkid1) {
      for(size_t shrid1=0; shrid1<shr_v.size(); ++shrid1) {
	flashana::QCluster_t qcluster_proton_shower  = build1e1pQCluster(trkid1, vtx, *(shr_v[shrid1]), dedxgen_v);
	qcluster_proton_shower_v.emplace_back(std::move(qcluster_proton_shower));
	
	pair_v[0] = trkid1;
	pair_v[1] = shrid1;
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
    // make flash hypothesis for proton-muon pair
    //
    for(const auto& qcluster_1mu1p : qcluster_proton_muon_v) {
      flashana::Flash_t hypo_1mu1p = genflashmatch->GenerateUnfittedFlashHypothesis(qcluster_1mu1p);

      proton_muon_chi2_1mu1p_v.resize(proton_muon_chi2_1mu1p_v.size()+1);
      auto& proton_muon_chi2_1mu1p = proton_muon_chi2_1mu1p_v.back();

      proton_muon_chi2_shape_1mu1p_v.resize(proton_muon_chi2_shape_1mu1p_v.size()+1);
      auto& proton_muon_chi2_shape_1mu1p = proton_muon_chi2_shape_1mu1p_v.back();

      proton_muon_hypo_totpe_v.resize(proton_muon_hypo_totpe_v.size()+1);
      auto& proton_muon_hypo_totpe = proton_muon_hypo_totpe_v.back();

      proton_muon_hypo_pe_vv.resize(proton_muon_hypo_pe_vv.size()+1);
      auto& proton_muon_hypo_pe_v = proton_muon_hypo_pe_vv.back();

      proton_muon_best_data_flash_v.resize(proton_muon_best_data_flash_v.size()+1);
      auto& proton_muon_best_data_flash = proton_muon_best_data_flash_v.back();

      proton_muon_hypo_pe_v.resize(32);

      FillChi2(dataflash_v,
	       hypo_1mu1p,
	       proton_muon_chi2_1mu1p,
	       proton_muon_chi2_shape_1mu1p,
	       proton_muon_hypo_totpe,
	       proton_muon_hypo_pe_v,
	       proton_muon_best_data_flash);
    }
    
    //
    // make flash hypothesis for muon-proton pair
    //
    for(const auto& qcluster_1mu1p : qcluster_muon_proton_v) {
      flashana::Flash_t hypo_1mu1p = genflashmatch->GenerateUnfittedFlashHypothesis(qcluster_1mu1p);
      
      muon_proton_chi2_1mu1p_v.resize(muon_proton_chi2_1mu1p_v.size()+1);
      auto& muon_proton_chi2_1mu1p = muon_proton_chi2_1mu1p_v.back();

      muon_proton_chi2_shape_1mu1p_v.resize(muon_proton_chi2_shape_1mu1p_v.size()+1);
      auto& muon_proton_chi2_shape_1mu1p = muon_proton_chi2_shape_1mu1p_v.back();

      muon_proton_hypo_totpe_v.resize(muon_proton_hypo_totpe_v.size()+1);
      auto& muon_proton_hypo_totpe = muon_proton_hypo_totpe_v.back();

      muon_proton_hypo_pe_vv.resize(muon_proton_hypo_pe_vv.size()+1);
      auto& muon_proton_hypo_pe_v = muon_proton_hypo_pe_vv.back();

      muon_proton_best_data_flash_v.resize(muon_proton_best_data_flash_v.size()+1);
      auto& muon_proton_best_data_flash = muon_proton_best_data_flash_v.back();

      muon_proton_hypo_pe_v.resize(32);

      FillChi2(dataflash_v,
	       hypo_1mu1p,
	       muon_proton_chi2_1mu1p,
	       muon_proton_chi2_shape_1mu1p,
	       muon_proton_hypo_totpe,
	       muon_proton_hypo_pe_v,
	       muon_proton_best_data_flash);
    }

    //
    // make flash hypothesis for proton-shower pair
    //
    for(const auto& qcluster_1e1p : qcluster_proton_shower_v) {
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

  void InterSelFlashMatch::FillChi2(const std::vector<flashana::Flash_t>& dataflash_v,
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

  flashana::QCluster_t InterSelFlashMatch::build1mu1pQCluster(const int protonid,
							      const int muonid,
							      std::vector<larlitecv::TrackHitSorter>& dedxgen_v) {
    // we use the proton id to build proton hypothesis
    // other tracks are assumed as muons

    const float ly_proton   = 19200;
    const float ly_muon     = 24000;
    
    flashana::QCluster_t qinteraction; // should reserve
    qinteraction.reserve(1000);
    
    for ( int itrack=0; itrack<(int)dedxgen_v.size(); itrack++ ) {

      float ly = kINVALID_FLOAT;
      if (itrack == protonid) 
	ly = ly_proton;
      else if (itrack == muonid) 
	ly = ly_muon;
      else 
	continue;

      // we get the dedx track
      std::vector<std::vector<float>> dedx_track_per_plane(3);
      dedxgen_v[itrack].getPathBinneddEdx(0.5,0.5,dedx_track_per_plane);

      // todo: use v plane if y plane missing too many pieces
      const std::vector< std::vector<float> >& bincenter_xyz = dedxgen_v[itrack].getBinCentersXYZ(2); 
      LLCV_DEBUG() << "1mu1p-track: bincenters=" << bincenter_xyz.size() << " vs " << dedx_track_per_plane[2].size() << std::endl;
      
      for (int ipt=0; ipt<(int)dedx_track_per_plane[2].size(); ipt++) {
	if ( ipt>=(int)bincenter_xyz.size() ) continue;
	
	float dedx = dedx_track_per_plane[2].at(ipt);

        if ( dedx<1.0e-2 )
 	  dedx = 2.07; // MeV/cm


	const std::vector<float>& edep_pos = bincenter_xyz[ipt];
	if ( edep_pos.size()==3 ) {
	  float numphotons = dedx*(2*0.5)*ly;
	  LLCV_DEBUG() << "1mu1p-track: (" << edep_pos[0] << "," << edep_pos[1] << "," << edep_pos[2] << ")" 
		       << " numphotons=" << numphotons << " dedx=" << dedx << std::endl;
	  flashana::QPoint_t qpt( edep_pos[0], edep_pos[1], edep_pos[2], numphotons );
	  qinteraction.emplace_back( std::move(qpt) );
	} // end edep_pos.size
      } // end dedx_track_per_plane[2]
    } // end track loop

    return qinteraction;
  }

  flashana::QCluster_t InterSelFlashMatch::build1e1pQCluster(const int protonid,
							     const larlite::vertex& vtx, 
							     const larlite::shower& shower, 
							     std::vector<larlitecv::TrackHitSorter>& dedxgen_v) {
    
    // we use the proton id to build proton hypothesis. same with shrid.
    // any leftover tracks are assumed as muons

    const float ly_proton = 19200;
    //const float ly_muon   = 24000;
    const float ly_shower = 20000;
    
    flashana::QCluster_t qinteraction; // should reserve
    qinteraction.reserve(1000);
    size_t itrack = protonid;
    
    // we get the dedx track
    std::vector< std::vector<float> > dedx_track_per_plane(3);
    dedxgen_v[itrack].getPathBinneddEdx( 0.5, 0.5, dedx_track_per_plane );
    // todo: use v plane if y plane missing too many pieces
    
    const std::vector< std::vector<float> >& bincenter_xyz = dedxgen_v[itrack].getBinCentersXYZ( 2 ); 
    LLCV_DEBUG() << "1e1p-track: bincenters=" << bincenter_xyz.size() << " vs " << dedx_track_per_plane[2].size() << std::endl;
    
    float ly = ly_proton;
    
    for (int ipt=0; ipt<(int)dedx_track_per_plane[2].size(); ipt++) {
      if ( ipt>=(int)bincenter_xyz.size() ) continue;
      float dedx = dedx_track_per_plane[2].at(ipt);
      const std::vector<float>& edep_pos = bincenter_xyz[ipt];
      if ( edep_pos.size()==3 ) {
	float numphotons = dedx*(2*0.5)*ly;
	flashana::QPoint_t qpt( edep_pos[0], edep_pos[1], edep_pos[2], numphotons );
	LLCV_DEBUG() << "1e1p-track: (" << edep_pos[0] << "," << edep_pos[1] << "," << edep_pos[2] << ") numphotons=" << numphotons << std::endl;
	qinteraction.emplace_back( std::move(qpt) );
      }
    }
    
    //make shower qcluser: each point corresponds to the number of photons
    
    float shower_energy = shower.Energy_v().at(2); 
    if (shower_energy<=0) {
      shower_energy = shower.Energy_v().at(0);
      shower_energy+= shower.Energy_v().at(1);
      shower_energy/=2.0;
    }
    float used_length = shower_energy/2.4; // MeV / (MeV/cm)
    
    LLCV_DEBUG() << " shower: shower energy=" << shower_energy << " used length=" << used_length << std::endl;
    
    float maxstepsize = 0.3;
    int nsteps = used_length/maxstepsize+1;
    float step = used_length/float(nsteps);
    for (int istep=0; istep<=nsteps; istep++) {
      double pos[3];
      pos[0] = vtx.X() + (step*istep)*shower.Direction().X();
      pos[1] = vtx.Y() + (step*istep)*shower.Direction().Y();
      pos[2] = vtx.Z() + (step*istep)*shower.Direction().Z();      
      float numphotons = (shower_correction_factor*step)*ly_shower;
      flashana::QPoint_t qpt( pos[0], pos[1], pos[2], numphotons );
      LLCV_DEBUG() << "1e1p-shower: (" << pos[0] << "," << pos[1] << "," << pos[2] << ") numphotons=" << numphotons << std::endl;
      qinteraction.push_back( qpt );
    }
    
    return qinteraction;
  }
  
  void InterSelFlashMatch::Finalize() {
    LLCV_DEBUG() << "start" << std::endl;

    _fout->cd();
    outtree->Write();

    LLCV_DEBUG() << "end" << std::endl;
  }

  void InterSelFlashMatch::ResetVertex() {
    
    vtxpos_x = kINVALID_FLOAT;
    vtxpos_y = kINVALID_FLOAT;
    vtxpos_z = kINVALID_FLOAT;

    number_tracks = -1*kINVALID_INT;
    number_showers= -1*kINVALID_INT;
    
    // data flash
    ndata_flashes = -1.0*kINVALID_INT;
    data_totpe_v.clear();
    data_pe_vv.clear();

    // track and track ID
    proton_muon_pair_id_vv.clear();
    proton_muon_chi2_1mu1p_v.clear();
    proton_muon_chi2_shape_1mu1p_v.clear();
    proton_muon_hypo_totpe_v.clear();
    proton_muon_hypo_pe_vv.clear();
    proton_muon_best_data_flash_v.clear();

    muon_proton_pair_id_vv.clear();
    muon_proton_chi2_1mu1p_v.clear();
    muon_proton_chi2_shape_1mu1p_v.clear();
    muon_proton_hypo_totpe_v.clear();
    muon_proton_hypo_pe_vv.clear();
    muon_proton_best_data_flash_v.clear();

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
