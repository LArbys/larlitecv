#ifndef __LARLITECV_LARBYSMC_CXX__
#define __LARLITECV_LARBYSMC_CXX__

#include "LArbysMC.h"

#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"
#include "Base/MCConstants.h"
#include "DataFormat/storage_manager.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctruth.h"

namespace larlitecv {
  
  //static LArbysImageMCProcessFactory __global_LArbysImageMCProcessFactory__;

  LArbysMC::LArbysMC()
    : _name("LArbysMC"),
      _mc_tree(nullptr)
  {
    _producer_mctruth  = "generator";
    _producer_mctrack  = "mcreco";
    _producer_mcshower = "mcreco";
    _psce = new larutil::SpaceChargeMicroBooNE();
  }

  LArbysMC::~LArbysMC()
  {
    if (_mc_tree)
      delete _mc_tree;
    _mc_tree = nullptr;
  }

  // void LArbysMC::configure(const PSet& cfg)
  // {
  //   // _rse_producer       = cfg.get<std::string>("RSEProducer");
  //   // _producer_roi       = cfg.get<std::string>("MCProducer");
  //   // _producer_image2d   = cfg.get<std::string>("Image2DProducer");
  //   // _neutrino_present   = cfg.get<bool>("NeutrinoPresent",true);
  // }

  /**
   * initialize interal storage tree
   *
   */
  void LArbysMC::initialize()
  {
    _mc_tree = new TTree("LArbysMCTree","MC infomation");
    bindAnaVariables(_mc_tree);
  }

  /**
   * bind variables to the TTree we will fill to
   *
   */
  void LArbysMC::bindAnaVariables( TTree* mc_tree )
  {

    // event indexing
    mc_tree->Branch("run",     &_run,     "run/I");
    mc_tree->Branch("subrun",  &_subrun,  "subrun/I");
    mc_tree->Branch("event",   &_event,   "event/I");
    mc_tree->Branch("entry"  , &_entry,   "entry/I");

    mc_tree->Branch("nuPresent",       &_neutrino_present, "nuPresent/O" );
    mc_tree->Branch("currentType",     &_current_type,     "currentType/I");
    mc_tree->Branch("interactionType", &_interaction_type, "interactionType/I");
    mc_tree->Branch("genieMode",       &_genie_mode,       "genieMode/I");
    mc_tree->Branch("Enu_true",        &_Enu_true,         "Enu_true/F");
    mc_tree->Branch("vtx_x",           &_vtx_x,            "vtx_x/F");
    mc_tree->Branch("vtx_y",           &_vtx_y,            "vtx_y/F");
    mc_tree->Branch("vtx_z",           &_vtx_z,            "vtx_z/F");        
    mc_tree->Branch("vtx_sce_x",       &_vtx_sce_x,        "vtx_sce_x/F");
    mc_tree->Branch("vtx_sce_y",       &_vtx_sce_y,        "vtx_sce_y/F");
    mc_tree->Branch("vtx_sce_z",       &_vtx_sce_z,        "vtx_sce_z/F");
    mc_tree->Branch("vtx_tick",        &_vtx_tick,         "vtx_tick/F");
    mc_tree->Branch("vtx_wire",        _vtx_wire,          "vtx_wire[3]/F");

    mc_tree->Branch("evis",            &_evis,             "evis/F");
    mc_tree->Branch("evis_had",        &_evis_had,         "evis_had/F");
    mc_tree->Branch("evis_vtx",        &_evis_vtx,         "evis_vtx/F");
    mc_tree->Branch("evis_lep",        &_evis_lep,         "evis_lep/F");

    mc_tree->Branch("hi_lep_pdg",      &_hi_lep_pdg,       "hi_lep_pdg/I");
    
    mc_tree->Branch("nprimary", &_nprimary,"nprimary/I");

    mc_tree->Branch("nproton", &_nproton,"nproton/I");
    mc_tree->Branch("nproton_60mev", &_nproton_60mev,"nproton_60mev/I");    
    mc_tree->Branch("nlepton", &_nlepton,"nlepton/I");
    mc_tree->Branch("nlepton_35mev", &_nlepton_35mev,"nlepton_35mev/I");    
    mc_tree->Branch("nmeson", &_nmeson,"nmeson/I");
    mc_tree->Branch("nmeson_35mev", &_nmeson_35mev,"nmeson_35mev/I");
    mc_tree->Branch("nshower", &_nshower,"nshower/I");
    mc_tree->Branch("nneutron", &_nneutron,"nneutron/I");
    mc_tree->Branch("npi0",     &_npi0, "npi0/I" );
    mc_tree->Branch("is1l1p0pi", &_1l1p0pi, "is1l1p0pi/I");
    
    
    // mc_tree->Branch("parentPDG",&_parent_pdg, "parentPDG/I");
    // mc_tree->Branch("signal",   &_is_signal,  "_is_signal/O");

    // mc_tree->Branch("energyDeposit",&_energy_deposit,"energyDeposit/D");
    // mc_tree->Branch("energyInit",&_energy_init,"energyInit/D");

    // mc_tree->Branch("parentX",&_parent_x,"parentX/D");
    // mc_tree->Branch("parentY",&_parent_y,"parentY/D");
    // mc_tree->Branch("parentZ",&_parent_z,"parentZ/D");
    // mc_tree->Branch("parentT",&_parent_t,"parentT/D");

    // mc_tree->Branch("parentPx",&_parent_px,"parentPx/D");
    // mc_tree->Branch("parentPy",&_parent_py,"parentPy/D");
    // mc_tree->Branch("parentPz",&_parent_pz,"parentpz/D");


    // mc_tree->Branch("vtx2d_w",&_vtx_2d_w_v);
    // mc_tree->Branch("vtx2d_t",&_vtx_2d_t_v);


    // mc_tree->Branch("hi_lep_pdg",&_hi_lep_pdg,"hi_lep_pdg/I");

    // mc_tree->Branch("dep_sum_lepton",&_dep_sum_lepton,"dep_sum_lepton/D");
    // mc_tree->Branch("ke_sum_lepton",&_ke_sum_lepton,"ke_sum_lepton/D");

    // mc_tree->Branch("dep_sum_proton",&_dep_sum_proton,"dep_sum_proton/D");
    // mc_tree->Branch("ke_sum_proton",&_ke_sum_proton,"ke_sum_proton/D");

    // mc_tree->Branch("dep_sum_meson",&_dep_sum_meson,"dep_sum_meson/D");
    // mc_tree->Branch("ke_sum_meson",&_ke_sum_meson,"ke_sum_meson/D");

    // mc_tree->Branch("dep_sum_shower",&_dep_sum_shower,"dep_sum_shower/D");
    // mc_tree->Branch("ke_sum_shower",&_ke_sum_shower,"ke_sum_shower/D");

    // mc_tree->Branch("dep_sum_neutron",&_dep_sum_neutron,"dep_sum_neutron/D");
    // mc_tree->Branch("ke_sum_neutron",&_ke_sum_neutron,"ke_sum_neutron/D");

    // mc_tree->Branch("daughterPx_v", &_daughterPx_v);
    // mc_tree->Branch("daughterPy_v", &_daughterPy_v);
    // mc_tree->Branch("daughterPz_v", &_daughterPz_v);

    // mc_tree->Branch("daughterX_v", &_daughterX_v);
    // mc_tree->Branch("daughterY_v", &_daughterY_v);
    // mc_tree->Branch("daughterZ_v", &_daughterZ_v);

    // mc_tree->Branch("daughterPdg_v"           , &_daughter_pdg_v);
    // mc_tree->Branch("daughterTrackid_v"       , &_daughter_trackid_v);
    // mc_tree->Branch("daughterParentTrackid_v" , &_daughter_parenttrackid_v);
    // mc_tree->Branch("daughterLength3d_v"      , &_daughter_length3d_v);
    // mc_tree->Branch("daughterLength_vv"       , &_daughter_length_vv);
    // mc_tree->Branch("daughterEnergyInit_v"    , &_daughter_energyinit_v);
    // mc_tree->Branch("daughterEnergyDep_v"     , &_daughter_energydep_v);
    // mc_tree->Branch("daughter2DStartX_vv"     , &_daughter_2dstartx_vv);
    // mc_tree->Branch("daughter2DStartY_vv"     , &_daughter_2dstarty_vv);
    // mc_tree->Branch("daughter2DEndX_vv"       , &_daughter_2dendx_vv);
    // mc_tree->Branch("daughter2DEndY_vv"       , &_daughter_2dendy_vv);
    // mc_tree->Branch("daughter2DCosAngle_vv"   , &_daughter_2dcosangle_vv);

    // mc_tree->Branch("true_proton_end_pt_x", &_true_proton_end_pt_x);
    // mc_tree->Branch("true_proton_end_pt_y", &_true_proton_end_pt_y);
    // mc_tree->Branch("true_proton_end_pt_z", &_true_proton_end_pt_z);

    // mc_tree->Branch("proton_1st_pt_x", &_proton_1st_pt_x);
    // mc_tree->Branch("proton_1st_pt_y", &_proton_1st_pt_y);
    // mc_tree->Branch("proton_1st_pt_z", &_proton_1st_pt_z);

    // mc_tree->Branch("proton_last_pt_x", &_proton_last_pt_x);
    // mc_tree->Branch("proton_last_pt_y", &_proton_last_pt_y);
    // mc_tree->Branch("proton_last_pt_z", &_proton_last_pt_z);

    // mc_tree->Branch("true_lepton_end_pt_x", &_true_lepton_end_pt_x);
    // mc_tree->Branch("true_lepton_end_pt_y", &_true_lepton_end_pt_y);
    // mc_tree->Branch("true_lepton_end_pt_z", &_true_lepton_end_pt_z);

    // mc_tree->Branch("lepton_1st_pt_x", &_lepton_1st_pt_x);
    // mc_tree->Branch("lepton_1st_pt_y", &_lepton_1st_pt_y);
    // mc_tree->Branch("lepton_1st_pt_z", &_lepton_1st_pt_z);

    // mc_tree->Branch("lepton_last_pt_x", &_lepton_last_pt_x);
    // mc_tree->Branch("lepton_last_pt_y", &_lepton_last_pt_y);
    // mc_tree->Branch("lepton_last_pt_z", &_lepton_last_pt_z);

  }

  bool LArbysMC::process(larlite::storage_manager& mgr)
  {
    
    Clear();

    larlite::event_mctruth* ev_mctruth =
      (larlite::event_mctruth*)mgr.get_data(larlite::data::kMCTruth,_producer_mctruth);

    larlite::event_mctrack* ev_mctrack =
      (larlite::event_mctrack*)mgr.get_data(larlite::data::kMCTrack,_producer_mctrack);

    larlite::event_mcshower* ev_mcshower =
      (larlite::event_mcshower*)mgr.get_data(larlite::data::kMCShower,_producer_mcshower);
    
    _run    = mgr.run_id();
    _subrun = mgr.subrun_id();
    _event  = mgr.event_id();
    _entry  = (int)mgr.get_index();

    // If we've got a neutrino, sample that
    auto const& mct = ev_mctruth->at(0);
    _neutrino_present = mct.Origin()==larlite::simb::kBeamNeutrino;
    
    // way to fill an empty entry?
    if ( !_neutrino_present ) {
      if ( _mc_tree )
        _mc_tree->Fill();
      return false;
    }
        
    auto const& nu = mct.GetNeutrino();

    _current_type     = nu.CCNC();
    _interaction_type = nu.InteractionType();
    _genie_mode       = nu.Mode();
    _Enu_true         = nu.Nu().Trajectory().front().E()*1000.0;
    _vtx_t            = nu.Nu().Trajectory().front().T();
    _vtx_x            = nu.Nu().Trajectory().front().X();
    _vtx_y            = nu.Nu().Trajectory().front().Y();
    _vtx_z            = nu.Nu().Trajectory().front().Z();    

    // SCE correction
    std::vector<double> pos_offset = _psce->GetPosOffsets( _vtx_x, _vtx_y, _vtx_z );
    _vtx_sce_x = _vtx_x - pos_offset[0] + 0.7;
    _vtx_sce_y = _vtx_y + pos_offset[1];
    _vtx_sce_z = _vtx_z + pos_offset[2];
    
    const float trig_time = 4050.0;
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    _vtx_tick = ( _vtx_t*1.0e-3 - (trig_time-4050.0) )/0.5 + _vtx_sce_x/cm_per_tick + 3200.0;

    if ( _vtx_y<-117.0 || _vtx_y>117.0
         || _vtx_z<0 || _vtx_z>1036.0 ) {
      // out of tpc
      for (int p=0; p<3; p++ ) _vtx_wire[p] = -1;
    }
    else {
      // in TPC
      double dpos[3] = { _vtx_sce_x, _vtx_sce_y, _vtx_sce_z };
      for (int p=0;p<3; p++) {
        _vtx_wire[p] = larutil::Geometry::GetME()->WireCoordinate( dpos, p );
      }
    }
        
    // _parent_pdg       = roi.PdgCode();
    // _energy_deposit   = roi.EnergyDeposit();
    // _energy_init      = roi.EnergyInit();
    // _parent_x         = roi.X();
    // _parent_y         = roi.Y();
    // _parent_z         = roi.Z();
    // _parent_t         = roi.T();
    // _parent_px        = roi.Px();
    // _parent_py        = roi.Py();
    // _parent_pz        = roi.Pz();

    // final state tally
    _nprimary = 0;
    _ntotal   = 0;
    _nproton  = 0;
    _nproton_60mev  = 0;    
    _nneutron = 0;
    _nlepton  = 0;
    _nlepton_35mev  = 0;    
    _nmeson   = 0;
    _nmeson_35mev   = 0;    
    _nshower  = 0;
    _npi0 = 0;
    _1l1p0pi = 0;
    
    _evis = 0;
    _evis_had = 0;
    _evis_vtx = 0;
    _evis_lep = 0;

    _hi_lep_pdg      = -1;
    _hi_lep_e        = 0;


    int vidx = -1;
    for(auto const& mct : *ev_mctrack ) {
      vidx++;
      if ( mct.Origin()==1 ) {
        // neutrino interaction primary
        int pid  = mct.PdgCode();

        std::cout << "track[" << vidx << "] origin=" << mct.Origin()
                  << " tid=" << mct.TrackID()
                  << " [ mid=" << mct.MotherTrackID()
                  << " mpdg=" << mct.MotherPdgCode() << " ]"
                  << " [ aid=" << mct.AncestorTrackID()
                  << " apdg=" << mct.AncestorPdgCode() << "]"
                  << " pid=" << mct.PdgCode()
                  << " E=" << mct.Start().E()          
                  << std::endl;
        
        
        if ( abs(pid)==12 || abs(pid)==14 || abs(pid)==16 )
          continue; // neutrino: skip

        if ( mct.TrackID()!=mct.MotherTrackID() )
          continue; // skip secondaries
        
        _nprimary++;
        
        float pE = mct.Start().E();
        if ( abs(pid)==13 ){
          // muon
          float ke=pE-105.0;
          if ( ke>35.0 ) _nlepton_35mev++;
          _nlepton++;
          _evis_lep += ke;
          _evis += ke;
          if ( ke>_hi_lep_e ) {
            _hi_lep_e = ke;
            _hi_lep_pdg =pid;
          }
        }        
        else if ( abs(pid)==11 ) {
          // electron
          float ke = pE-0.511;
          if ( ke>35.0 ) _nlepton_35mev++;
          _nlepton++;
          _evis += ke;
          _evis_lep += ke;
          if ( ke>_hi_lep_e ) {
            _hi_lep_e = ke;
            _hi_lep_pdg =pid;
          }          
        }
        else if ( pid==2212 ) {
          // proton
          float ke = pE-938.0;
          if( ke>60.0 ) {
            _nproton_60mev++;
          }
          else {
            _evis_vtx += ke;
          }
          _nproton++;
          _evis += ke;
          _evis_had += ke;
        }
        else if ( abs(pid)==211 ) {
          float ke = pE-135.0;
          if ( ke>35.0 ) {
            _nmeson_35mev++;
          }
          else {
            _evis_vtx += ke;
          }
          _nmeson++;          
          _evis_had += ke;
          _evis += ke;
        }
        else if ( pid==2112 ) {
          _nneutron++;
        }
        
      }
    }

    vidx = -1;
    std::set<int> pi0_ids; // list of pi0 id's we've already used
    for(auto const& mcs : *ev_mcshower ) {
      vidx++;
      if ( mcs.Origin()==1 ) {
        // neutrino interaction primary
        int pid  = mcs.PdgCode();

        std::cout << "shower[" << vidx << "] origin=" << mcs.Origin()
                  << " tid=" << mcs.TrackID()
                  << " [ mid=" << mcs.MotherTrackID()
                  << " mpd=" << mcs.MotherPdgCode() << " ]"
                  << " [ aid=" << mcs.AncestorTrackID()
                  << " apdg=" << mcs.AncestorPdgCode() << " ]"
                  << " pid=" << mcs.PdgCode()
                  << " E=" << mcs.Start().E()
                  << std::endl;
        
        if ( abs(pid)==12 || abs(pid)==14 || abs(pid)==16 )
          continue; // neutrino: skip

        if ( abs(pid)!=22 && mcs.TrackID()!=mcs.MotherTrackID() )
          continue; // skip secondaries (non photons)
        if ( abs(pid)==22 && mcs.MotherTrackID()!=mcs.AncestorTrackID() )
          continue; // skip secondaries (photons)
        
        float pE = mcs.Start().E();
        if ( abs(pid)==11 ) {
          // electron
          float ke = pE-0.511;
          if ( ke>35.0 ) _nlepton_35mev++;
          _nlepton++;
          _nshower++;          
          _evis += ke;
          _evis_lep += ke;
          if ( ke>_hi_lep_e ) {
            _hi_lep_e = ke;
            _hi_lep_pdg =pid;
          }          
        }
        else if ( abs(pid)==22 ) {

          if ( mcs.MotherPdgCode()==111 ) {
            // this is a pi0 photon
            int mid = mcs.MotherTrackID();
            if (  pi0_ids.find(mid)==pi0_ids.end() ) {
              // new id
              pi0_ids.insert(mid);
              _npi0++;
              _nmeson++;
            }
          }
          
          // photon
          float ke = pE;
          _nshower++;
          _evis += ke;
        }
        
      }
    }

    if ( _nproton_60mev==1 && _nlepton_35mev==1 && _nmeson_35mev==0 && _npi0==0 )
      _1l1p0pi = 1;
    
    // _dep_sum_lepton  = 0;
    // _ke_sum_lepton   = 0;

    // _dep_sum_proton  = 0;
    // _ke_sum_proton   = 0;

    // _dep_sum_meson   = 0;
    // _ke_sum_meson    = 0;

    // _dep_sum_shower  = 0;
    // _ke_sum_shower   = 0;

    // _dep_sum_neutron = 0;
    // _ke_sum_neutron  = 0;

  
    // bool first=_neutrino_present;
    // std::vector<aparticle> proton_v;
    // proton_v.clear();
    
    // // loop through mctrack final state, looking for neutrino origin only
    // for(auto const& mct : *ev_mctrack ) {

    //   if (first){
    //     // assumes the first entry is the neutrino ...
    //     first=false;
    //     continue;
    //   }

      
    //   std::cout << "This particle is PDG code " << roi.ParentPdgCode() << std::endl;

    //   // Get a unit vector for this pdg in 3 coordinates
    //   auto px = roi.Px();
    //   auto py = roi.Py();
    //   auto pz = roi.Pz();

    //   auto lenp = sqrt(px*px+py*py+pz*pz);

    //   px /= lenp;
    //   py /= lenp;
    //   pz /= lenp;

    //   // The Original location
    //   auto x0 = roi.X();
    //   auto y0 = roi.Y();
    //   auto z0 = roi.Z();
    //   auto t0 = roi.T();

    //   // Here is another point in the direction of p, ignore units...
    //   // Pxyz are info from genie (meaning that it won't be identical to PCA assumption).
    //   auto end_pt = roi.EndPosition();
    //   auto first_step = roi.FirstStep();
    //   auto last_step = roi.LastStep();      

    //   auto x1 = x0 + px;
    //   auto y1 = y0 + py;
    //   auto z1 = z0 + pz;
    //   auto t1 = t0;

    //   _daughter_length_vv.resize(3);
    //   _daughter_2dstartx_vv.resize(3);
    //   _daughter_2dstarty_vv.resize(3);
    //   _daughter_2dendx_vv.resize(3);
    //   _daughter_2dendy_vv.resize(3);
    //   _daughter_2dcosangle_vv.resize(3);

    //   //lets project both points

    //   for(size_t plane=0; plane<3; ++plane) {

    //     auto& daughter_length_v      = _daughter_length_vv[plane];
    //     auto& daughter_2dstartx_v    = _daughter_2dstartx_vv[plane];
    //     auto& daughter_2dstarty_v    = _daughter_2dstarty_vv[plane];
    //     auto& daughter_2dendx_v      = _daughter_2dendx_vv[plane];
    //     auto& daughter_2dendy_v      = _daughter_2dendy_vv[plane];
    //     auto& daughter_2dcosangle_v  = _daughter_2dcosangle_vv[plane];

    //     const auto& img  = ev_image2d->Image2DArray()[plane];
    //     const auto& meta = img.meta();

    //     double x_pixel0, y_pixel0, x_pixel1, y_pixel1;
    //     x_pixel0 = y_pixel0 = 0;
    //     x_pixel1 = y_pixel1 = 0;

    //     Project3D(meta,x0,y0,z0,t0,plane,x_pixel0,y_pixel0);
    //     Project3D(meta,x1,y1,z1,t1,plane,x_pixel1,y_pixel1);

    //     // Start and end in 2D
    //     geo2d::Vector<float> start(x_pixel0,y_pixel0);
    //     geo2d::Vector<float> end  (x_pixel1,y_pixel1);

    //     // Get the line of particle
    //     geo2d::HalfLine<float> hline(start,end-start);

    //     // Here is the bbox on this plane --> we need to get the single intersection point
    //     // for the half line and this bbox
    //     ImageMeta bb;
    //     try { bb = roi.BB(plane); }
    //     catch(...) // No intersection occurs?
    //       { continue; }

    //     auto roi_on_plane = Get2DRoi(meta,roi.BB(plane));

    //     // The start point will be inside the 2D ROI
    //     // we need to intersection point between the edge and this half line, find it

    //     auto pt_edge = Intersection(hline,roi_on_plane);

    //     daughter_length_v.push_back(geo2d::length(pt_edge - start));

    //     daughter_2dstartx_v.push_back(start.x);
    //     daughter_2dstarty_v.push_back(start.y);
    //     daughter_2dendx_v.push_back(pt_edge.x);
    //     daughter_2dendy_v.push_back(pt_edge.y);

    //     auto dir = pt_edge - start;
    //     double cosangle = dir.x / sqrt(dir.x*dir.x + dir.y*dir.y);

    //     daughter_2dcosangle_v.push_back(cosangle);
    //   }
      
    //   _daughterPx_v.push_back(roi.Px());
    //   _daughterPy_v.push_back(roi.Py());
    //   _daughterPz_v.push_back(roi.Pz());
    //   _daughterX_v.push_back(roi.X());
    //   _daughterY_v.push_back(roi.Y());
    //   _daughterZ_v.push_back(roi.Z());
      
    //   double length3d = sqrt(pow(roi.LastStep().X()-roi.FirstStep().X(),2)+pow(roi.LastStep().Y()-roi.FirstStep().Y(),2)+pow(roi.LastStep().Z()-roi.FirstStep().Z(),2));
    //   _daughter_length3d_v.push_back(length3d);
      
    //   int pdgcode = roi.PdgCode();

    //   _daughter_pdg_v.push_back((int) roi.PdgCode());
    //   _daughter_trackid_v.push_back((int) roi.TrackID());
    //   _daughter_parenttrackid_v.push_back((int) roi.ParentTrackID());
    //   _daughter_energyinit_v.push_back(roi.EnergyInit());
    //   _daughter_energydep_v.push_back(roi.EnergyDeposit());

    //   _ntotal++;

    //   pdgcode = std::abs(pdgcode);

    //   // This is proton
    //   if (pdgcode==2212) {
    //     //primary protons
    //     if (roi.TrackID() == roi.ParentTrackID()) {
    //       _nproton++;
    //       _ke_sum_proton  += (roi.EnergyInit() - 938.0);
    //     }

    //     _true_proton_end_pt_x = end_pt.X();
    //     _true_proton_end_pt_y = end_pt.Y();
    //     _true_proton_end_pt_z = end_pt.Z();
	
    //     _proton_1st_pt_x = first_step.X();
    //     _proton_1st_pt_y = first_step.Y();
    //     _proton_1st_pt_z = first_step.Z();
	
    //     _proton_last_pt_x = last_step.X();
    //     _proton_last_pt_y = last_step.Y();
    //     _proton_last_pt_z = last_step.Z();      
	
    //     //capture proton
    //     aparticle thispro;
    //     thispro.trackid       = roi.TrackID();
    //     thispro.parenttrackid = roi.ParentTrackID();
    //     thispro.depeng        = roi.EnergyDeposit();
    //     thispro.primary = (roi.TrackID()==roi.ParentTrackID());
    //     proton_v.emplace_back(std::move(thispro));
    //   }
    //   // its not a primary, skip
    //   if (roi.TrackID() != roi.ParentTrackID()) continue;

    //   //it is
    //   _nprimary++;

    //   //this is neutron
    //   if (pdgcode==2112) {
    //     _nneutron++;
    //     _dep_sum_neutron += roi.EnergyDeposit();
    //     _ke_sum_neutron  += (roi.EnergyInit() - 939.5);
    //   }

    //   // mesons are pion,kaon, ... what else
    //   if (pdgcode==211 or pdgcode==321) {
    //     _nmeson++;
    //     _dep_sum_meson += roi.EnergyDeposit();
    //     _ke_sum_meson  += roi.EnergyInit();
    //   }

    //   // leptons are electron, muon, ... what else
    //   if (pdgcode==11 or pdgcode==13) {
    //     _nlepton++;
    //     _dep_sum_lepton += roi.EnergyDeposit();
    //     _ke_sum_lepton  += roi.EnergyInit();
    //     if (roi.EnergyInit() > _hi_lep_e)  {
    //       _hi_lep_e = roi.EnergyInit();
    //       _hi_lep_pdg = pdgcode;
    //     }
	
    //     _true_lepton_end_pt_x = end_pt.X();
    //     _true_lepton_end_pt_y = end_pt.Y();
    //     _true_lepton_end_pt_z = end_pt.Z();
	
    //     _lepton_1st_pt_x = first_step.X();
    //     _lepton_1st_pt_y = first_step.Y();
    //     _lepton_1st_pt_z = first_step.Z();
	
    //     _lepton_last_pt_x = last_step.X();
    //     _lepton_last_pt_y = last_step.Y();
    //     _lepton_last_pt_z = last_step.Z();      
	
    //   }

    //   // shower are electron, gamma, pi0
    //   if (pdgcode==11 or pdgcode==22 or pdgcode==111) {
    //     _nshower++;
    //     _dep_sum_shower += roi.EnergyDeposit();
    //     _ke_sum_shower  += roi.EnergyInit();
    //   }
    // }

    // auto nprotons=proton_v.size();

    // float max_proton_e = 0;
    // std::vector<float> proton_engs;
    // proton_engs.reserve(nprotons);
    

    // if (nprotons){

    //   for (int x=0; x < (int)nprotons; x++){
    //     const auto& proton1 = proton_v[x];
    //     max_proton_e  = 0;
    //     max_proton_e += proton1.depeng;
    //     auto trackid  = proton1.trackid;

    //     for (int y=x+1; y < (int)nprotons; y++ ){
    //       const auto& proton2 = proton_v[y];
    //       if (proton2.daughterof(proton1)) {
    //         max_proton_e += proton2.depeng;
    //       }
    //     }

    //     proton_engs.push_back(max_proton_e);
    //   }

    //   auto max_proton_e = *(std::max_element(std::begin(proton_engs),
    //     				     std::end(proton_engs)));
    //   _dep_sum_proton = max_proton_e;
    // }

    if ( _mc_tree )
      _mc_tree->Fill();

    return true;
  }


  void LArbysMC::finalize()
  {
    if ( _mc_tree )
      _mc_tree->Write();
  }

  void LArbysMC::Clear() {

    _neutrino_present = false;
    _current_type = -1;
    _interaction_type = -1;
    _genie_mode = -1;
    
    _vtx_2d_w_v.clear();
    _vtx_2d_t_v.clear();
    _vtx_2d_w_v.resize(3);
    _vtx_2d_t_v.resize(3);

    // _image_v.clear();
    // _image_v.resize(3);

    _daughter_pdg_v.clear();
    _daughter_trackid_v.clear();
    _daughter_parenttrackid_v.clear();

    _daughter_energyinit_v.clear();
    _daughter_energydep_v.clear();

    _daughter_length_vv.clear();
    _daughter_length3d_v.clear();

    _daughter_2dstartx_vv.clear();
    _daughter_2dstarty_vv.clear();

    _daughter_2dendx_vv.clear();
    _daughter_2dendy_vv.clear();

    _daughterPx_v.clear();
    _daughterPy_v.clear();
    _daughterPz_v.clear();

    _daughterX_v.clear();
    _daughterY_v.clear();
    _daughterZ_v.clear();

    _daughter_2dcosangle_vv.clear();
  }


  void LArbysMC::printInteractionInfo() {
    std::cout << "[LArbysMC::printInteractionInfo] ==============================================" << std::endl;
    std::cout << " rse: (" << _run << "," << _subrun << "," << _event << ") entry=" << _entry << std::endl;
    std::cout << " True Nu E: " << _Enu_true << " GeV" << std::endl;
    std::cout << " genie mode: " << _genie_mode << std::endl;
    std::cout << " is CC: " << (_current_type==0) << std::endl;
    std::cout << " interaction type: " << _interaction_type << std::endl;
    std::cout << " nlepton: " << _nlepton << "; above 35 MeV " << _nlepton_35mev << std::endl;
    std::cout << " nproton: " << _nproton << "; above 60 MeV " << _nproton_60mev << std::endl;
    std::cout << " nmeson: "  << _nmeson  << "; above 35 MeV " << _nmeson_35mev << std::endl;
    std::cout << " nneutron: " << _nneutron << std::endl;
    std::cout << " npi0: " << _npi0 << std::endl;
    std::cout << " is 1L1P (+ 0 pi + 0 pi0): " << _1l1p0pi << std::endl;
    std::cout << " evis: " << _evis << "; lepton: " << _evis_lep << "; hadronic: " << _evis_had << "; vertex " << _evis_vtx << std::endl;
    std::cout << " (x,y,z) true: (" << _vtx_x << "," << _vtx_y << "," << _vtx_z << ")" << std::endl;
    std::cout << " (x,y,z) sce: (" << _vtx_sce_x << "," << _vtx_sce_y << "," << _vtx_sce_z << ")" << std::endl;
    std::cout << " tick: " << _vtx_tick << std::endl;
    std::cout << " wires: (" << _vtx_wire[0] << "," << _vtx_wire[1] << "," << _vtx_wire[2] << ")" << std::endl;
    std::cout << "==============================================================================" << std::endl;    
  }
}
#endif
