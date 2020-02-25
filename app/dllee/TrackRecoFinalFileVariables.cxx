#include "TrackRecoFinalFileVariables.h"

namespace larlitecv {
namespace dllee {

  TrackRecoFinalFileVariables::TrackRecoFinalFileVariables( TTree* input_tree, TTree* output_tree ) {

    clearVariables();
    defineInputTreeVariables( input_tree );
    defineOutputTreeVariables( output_tree );    
    
  }

  TrackRecoFinalFileVariables::~TrackRecoFinalFileVariables() {
  }

  void TrackRecoFinalFileVariables::clearVariables() {

    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    trk_id_v = 0;
    trackQ3_v = 0;
    trackQ5_v = 0;
    trackQ10_v = 0;
    trackQ20_v = 0;
    trackQ30_v = 0;
    trackQ50_v = 0;
    E_muon_v = 0;
    E_proton_v = 0;
    Length_v = 0;
    Avg_Ion_v = 0;
    Avg_IonY_v = 0;
    Ion_5cm_v = 0;
    Ion_10cm_v = 0;
    Ion_tot_v = 0;
    IonY_5cm_v = 0;
    IonY_10cm_v = 0;
    IonY_tot_v = 0;
    Truncated_dQdX1_v = 0;
    Truncated_dQdX3_v = 0;
    IondivLength_v = 0;
    TotalADCvalues_v = 0;
    Angle_v = 0;
    Reco_goodness_v = 0;
    track_Goodness_v = 0;
    RecoVertex = 0;
    vertexPhi = 0;
    vertexTheta = 0;
    vertexPhi_2cm = 0;
    vertexTheta_2cm = 0;
    vertexPhi_5cm = 0;
    vertexTheta_5cm = 0;
    vertexPhi_7cm = 0;
    vertexTheta_7cm = 0;
    vertexPhi_10cm = 0;
    vertexTheta_10cm = 0;
    vertexPhi_12cm = 0;
    vertexTheta_12cm = 0;
    vertexPhi_15cm = 0;
    vertexTheta_15cm = 0;
    vertexPhi_17cm = 0;
    vertexTheta_17cm = 0;
    vertexPhi_20cm = 0;
    vertexTheta_20cm = 0;
    vertexPhi_30cm = 0;
    vertexTheta_30cm = 0;
    trackAvg15cm_x_sceadded = 0;
    trackAvg15cm_y_sceadded = 0;
    trackAvg15cm_z_sceadded = 0;
    closestWall = 0;
    recoEndPoints_x = 0;
    recoEndPoints_y = 0;
    recoEndPoints_z = 0;
    MCvertex = 0;
    
  }

  void TrackRecoFinalFileVariables::defineInputTreeVariables( TTree* input_tree ) {

    // need tracker trees
    input_tree->SetBranchAddress("run", &run, &b__run);
    input_tree->SetBranchAddress("subrun", &subrun, &b__subrun);
    input_tree->SetBranchAddress("event", &event, &b__event);
    input_tree->SetBranchAddress("entry", &entry, &b__entry);
    input_tree->SetBranchAddress("vtx_id", &vtx_id, &b_vtx_id);
    input_tree->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x_sceadded);
    input_tree->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y_sceadded);
    input_tree->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z_sceadded);
    input_tree->SetBranchAddress("randomSeed", &randomSeed, &b_randomSeed);
    input_tree->SetBranchAddress("trk_id_v", &trk_id_v, &b_trk_id_v);
    input_tree->SetBranchAddress("NtracksReco", &NtracksReco, &b_NtracksReco);
    input_tree->SetBranchAddress("MinLengthAllowed", &MinLengthAllowed, &b_MinLengthAllowed);
    input_tree->SetBranchAddress("trackQ3_v", &trackQ3_v, &b_trackQ3_v);
    input_tree->SetBranchAddress("trackQ5_v", &trackQ5_v, &b_trackQ5_v);
    input_tree->SetBranchAddress("trackQ10_v", &trackQ10_v, &b_trackQ10_v);
    input_tree->SetBranchAddress("trackQ20_v", &trackQ20_v, &b_trackQ20_v);
    input_tree->SetBranchAddress("trackQ30_v", &trackQ30_v, &b_trackQ30_v);
    input_tree->SetBranchAddress("trackQ50_v", &trackQ50_v, &b_trackQ50_v);
    input_tree->SetBranchAddress("E_muon_v", &E_muon_v, &b_E_muon_v);
    input_tree->SetBranchAddress("E_proton_v", &E_proton_v, &b_E_proton_v);
    input_tree->SetBranchAddress("Length_v", &Length_v, &b_Length_v);
    input_tree->SetBranchAddress("Avg_Ion_v", &Avg_Ion_v, &b_Avg_Ion_v);
    input_tree->SetBranchAddress("Avg_IonY_v", &Avg_IonY_v, &b_Avg_IonY_v);
    input_tree->SetBranchAddress("Ion_5cm_v", &Ion_5cm_v, &b_Ion_5cm_v);
    input_tree->SetBranchAddress("Ion_10cm_v", &Ion_10cm_v, &b_Ion_10cm_v);
    input_tree->SetBranchAddress("Ion_tot_v", &Ion_tot_v, &b_Ion_tot_v);
    input_tree->SetBranchAddress("IonY_5cm_v", &IonY_5cm_v, &b_IonY_5cm_v);
    input_tree->SetBranchAddress("IonY_10cm_v", &IonY_10cm_v, &b_IonY_10cm_v);
    input_tree->SetBranchAddress("IonY_tot_v", &IonY_tot_v, &b_IonY_tot_v);
    input_tree->SetBranchAddress("Truncated_dQdX1_v", &Truncated_dQdX1_v, &b_Truncated_dQdX1_v);
    input_tree->SetBranchAddress("Truncated_dQdX3_v", &Truncated_dQdX3_v, &b_Truncated_dQdX3_v);
    input_tree->SetBranchAddress("IondivLength_v", &IondivLength_v, &b_IondivLength_v);
    input_tree->SetBranchAddress("TotalADCvalues_v", &TotalADCvalues_v, &b_TotalADCvalues_v);
    input_tree->SetBranchAddress("Angle_v", &Angle_v, &b_Angle_v);
    input_tree->SetBranchAddress("Reco_goodness_v", &Reco_goodness_v, &b_Reco_goodness_v);
    input_tree->SetBranchAddress("track_Goodness_v", &track_Goodness_v, &b_track_Goodness_v);
    input_tree->SetBranchAddress("GoodVertex", &GoodVertex, &b_GoodVertex);
    input_tree->SetBranchAddress("missingTrack", &missingTrack, &b_missingTrack);
    input_tree->SetBranchAddress("nothingReconstructed", &nothingReconstructed, &b_nothingReconstructed);
    input_tree->SetBranchAddress("tooShortDeadWire", &tooShortDeadWire, &b_tooShortDeadWire);
    input_tree->SetBranchAddress("tooShortFaintTrack", &tooShortFaintTrack, &b_tooShortFaintTrack);
    input_tree->SetBranchAddress("tooManyTracksAtVertex", &tooManyTracksAtVertex, &b_tooManyTracksAtVertex);
    input_tree->SetBranchAddress("possibleCosmic", &possibleCosmic, &b_possibleCosmic);
    input_tree->SetBranchAddress("possiblyCrossing", &possiblyCrossing, &b_possiblyCrossing);
    input_tree->SetBranchAddress("branchingTracks", &branchingTracks, &b_branchingTracks);
    input_tree->SetBranchAddress("jumpingTracks", &jumpingTracks, &b_jumpingTracks);
    input_tree->SetBranchAddress("RecoVertex", &RecoVertex, &b_RecoVertex);
    input_tree->SetBranchAddress("vertexPhi", &vertexPhi, &b_vertexPhi);
    input_tree->SetBranchAddress("vertexTheta", &vertexTheta, &b_vertexTheta);
    input_tree->SetBranchAddress("vertexPhi_2cm", &vertexPhi_2cm, &b_vertexPhi_2cm);
    input_tree->SetBranchAddress("vertexTheta_2cm", &vertexTheta_2cm, &b_vertexTheta_2cm);
    input_tree->SetBranchAddress("vertexPhi_5cm", &vertexPhi_5cm, &b_vertexPhi_5cm);
    input_tree->SetBranchAddress("vertexTheta_5cm", &vertexTheta_5cm, &b_vertexTheta_5cm);
    input_tree->SetBranchAddress("vertexPhi_7cm", &vertexPhi_7cm, &b_vertexPhi_7cm);
    input_tree->SetBranchAddress("vertexTheta_7cm", &vertexTheta_7cm, &b_vertexTheta_7cm);
    input_tree->SetBranchAddress("vertexPhi_10cm", &vertexPhi_10cm, &b_vertexPhi_10cm);
    input_tree->SetBranchAddress("vertexTheta_10cm", &vertexTheta_10cm, &b_vertexTheta_10cm);
    input_tree->SetBranchAddress("vertexPhi_12cm", &vertexPhi_12cm, &b_vertexPhi_12cm);
    input_tree->SetBranchAddress("vertexTheta_12cm", &vertexTheta_12cm, &b_vertexTheta_12cm);
    input_tree->SetBranchAddress("vertexPhi_15cm", &vertexPhi_15cm, &b_vertexPhi_15cm);
    input_tree->SetBranchAddress("vertexTheta_15cm", &vertexTheta_15cm, &b_vertexTheta_15cm);
    input_tree->SetBranchAddress("vertexPhi_17cm", &vertexPhi_17cm, &b_vertexPhi_17cm);
    input_tree->SetBranchAddress("vertexTheta_17cm", &vertexTheta_17cm, &b_vertexTheta_17cm);
    input_tree->SetBranchAddress("vertexPhi_20cm", &vertexPhi_20cm, &b_vertexPhi_20cm);
    input_tree->SetBranchAddress("vertexTheta_20cm", &vertexTheta_20cm, &b_vertexTheta_20cm);
    input_tree->SetBranchAddress("vertexPhi_30cm", &vertexPhi_30cm, &b_vertexPhi_30cm);
    input_tree->SetBranchAddress("vertexTheta_30cm", &vertexTheta_30cm, &b_vertexTheta_30cm);
    input_tree->SetBranchAddress("trackAvg15cm_x_sceadded", &trackAvg15cm_x_sceadded, &b_trackAvg15cm_x_sceadded);
    input_tree->SetBranchAddress("trackAvg15cm_y_sceadded", &trackAvg15cm_y_sceadded, &b_trackAvg15cm_y_sceadded);
    input_tree->SetBranchAddress("trackAvg15cm_z_sceadded", &trackAvg15cm_z_sceadded, &b_trackAvg15cm_z_sceadded);
    input_tree->SetBranchAddress("closestWall", &closestWall, &b_closestWall);
    input_tree->SetBranchAddress("recoEndPoints_x", &recoEndPoints_x, &b_recoEndPoints_x);
    input_tree->SetBranchAddress("recoEndPoints_y", &recoEndPoints_y, &b_recoEndPoints_y);
    input_tree->SetBranchAddress("recoEndPoints_z", &recoEndPoints_z, &b_recoEndPoints_z);
    input_tree->SetBranchAddress("MCvertex", &MCvertex, &b_MCvertex);
    input_tree->SetBranchAddress("Ep_t", &Ep_t, &b_Ep_t);
    input_tree->SetBranchAddress("Em_t", &Em_t, &b_Em_t);
    input_tree->SetBranchAddress("Ee_t", &Ee_t, &b_Ee_t);
    input_tree->SetBranchAddress("MuonStartPoint_X", &MuonStartPoint_X, &b_MuonStartPoint_X);
    input_tree->SetBranchAddress("ProtonStartPoint_X", &ProtonStartPoint_X, &b_ProtonStartPoint_X);
    input_tree->SetBranchAddress("ElectronStartPoint_X", &ElectronStartPoint_X, &b_ElectronStartPoint_X);
    input_tree->SetBranchAddress("MuonStartPoint_Y", &MuonStartPoint_Y, &b_MuonStartPoint_Y);
    input_tree->SetBranchAddress("ProtonStartPoint_Y", &ProtonStartPoint_Y, &b_ProtonStartPoint_Y);
    input_tree->SetBranchAddress("ElectronStartPoint_Y", &ElectronStartPoint_Y, &b_ElectronStartPoint_Y);
    input_tree->SetBranchAddress("MuonStartPoint_Z", &MuonStartPoint_Z, &b_MuonStartPoint_Z);
    input_tree->SetBranchAddress("ProtonStartPoint_Z", &ProtonStartPoint_Z, &b_ProtonStartPoint_Z);
    input_tree->SetBranchAddress("ElectronStartPoint_Z", &ElectronStartPoint_Z, &b_ElectronStartPoint_Z);
    input_tree->SetBranchAddress("MuonEndPoint_X", &MuonEndPoint_X, &b_MuonEndPoint_X);
    input_tree->SetBranchAddress("ProtonEndPoint_X", &ProtonEndPoint_X, &b_ProtonEndPoint_X);
    input_tree->SetBranchAddress("ElectronEndPoint_X", &ElectronEndPoint_X, &b_ElectronEndPoint_X);
    input_tree->SetBranchAddress("MuonEndPoint_Y", &MuonEndPoint_Y, &b_MuonEndPoint_Y);
    input_tree->SetBranchAddress("ProtonEndPoint_Y", &ProtonEndPoint_Y, &b_ProtonEndPoint_Y);
    input_tree->SetBranchAddress("ElectronEndPoint_Y", &ElectronEndPoint_Y, &b_ElectronEndPoint_Y);
    input_tree->SetBranchAddress("MuonEndPoint_Z", &MuonEndPoint_Z, &b_MuonEndPoint_Z);
    input_tree->SetBranchAddress("ProtonEndPoint_Z", &ProtonEndPoint_Z, &b_ProtonEndPoint_Z);
    input_tree->SetBranchAddress("ElectronEndPoint_Z", &ElectronEndPoint_Z, &b_ElectronEndPoint_Z);
      
  }
  
  void TrackRecoFinalFileVariables::defineOutputTreeVariables( TTree* output_tree ) {

    // need tracker trees
    output_tree->Branch("run", &run);
    output_tree->Branch("subrun", &subrun);
    output_tree->Branch("event", &event);
    output_tree->Branch("entry", &entry);
    output_tree->Branch("vtx_id", &vtx_id);
    output_tree->Branch("vtx_x", &vtx_x);
    output_tree->Branch("vtx_y", &vtx_y);
    output_tree->Branch("vtx_z", &vtx_z);
    output_tree->Branch("randomSeed", &randomSeed);
    output_tree->Branch("trk_id_v", &trk_id_v);
    output_tree->Branch("NtracksReco", &NtracksReco);
    output_tree->Branch("MinLengthAllowed", &MinLengthAllowed);
    output_tree->Branch("trackQ3_v", &trackQ3_v);
    output_tree->Branch("trackQ5_v", &trackQ5_v);
    output_tree->Branch("trackQ10_v", &trackQ10_v);
    output_tree->Branch("trackQ20_v", &trackQ20_v);
    output_tree->Branch("trackQ30_v", &trackQ30_v);
    output_tree->Branch("trackQ50_v", &trackQ50_v);
    output_tree->Branch("E_muon_v", &E_muon_v);
    output_tree->Branch("E_proton_v", &E_proton_v);
    output_tree->Branch("Length_v", &Length_v);
    output_tree->Branch("Avg_Ion_v", &Avg_Ion_v);
    output_tree->Branch("Avg_IonY_v", &Avg_IonY_v);
    output_tree->Branch("Ion_5cm_v", &Ion_5cm_v);
    output_tree->Branch("Ion_10cm_v", &Ion_10cm_v);
    output_tree->Branch("Ion_tot_v", &Ion_tot_v);
    output_tree->Branch("IonY_5cm_v", &IonY_5cm_v);
    output_tree->Branch("IonY_10cm_v", &IonY_10cm_v);
    output_tree->Branch("IonY_tot_v", &IonY_tot_v);
    output_tree->Branch("Truncated_dQdX1_v", &Truncated_dQdX1_v);
    output_tree->Branch("Truncated_dQdX3_v", &Truncated_dQdX3_v);
    output_tree->Branch("IondivLength_v", &IondivLength_v);
    output_tree->Branch("TotalADCvalues_v", &TotalADCvalues_v);
    output_tree->Branch("Angle_v", &Angle_v);
    output_tree->Branch("Reco_goodness_v", &Reco_goodness_v);
    output_tree->Branch("track_Goodness_v", &track_Goodness_v);
    output_tree->Branch("GoodVertex", &GoodVertex);
    output_tree->Branch("missingTrack", &missingTrack);
    output_tree->Branch("nothingReconstructed", &nothingReconstructed);
    output_tree->Branch("tooShortDeadWire", &tooShortDeadWire);
    output_tree->Branch("tooShortFaintTrack", &tooShortFaintTrack);
    output_tree->Branch("tooManyTracksAtVertex", &tooManyTracksAtVertex);
    output_tree->Branch("possibleCosmic", &possibleCosmic);
    output_tree->Branch("possiblyCrossing", &possiblyCrossing);
    output_tree->Branch("branchingTracks", &branchingTracks);
    output_tree->Branch("jumpingTracks", &jumpingTracks);
    output_tree->Branch("RecoVertex", &RecoVertex);
    output_tree->Branch("vertexPhi", &vertexPhi);
    output_tree->Branch("vertexTheta", &vertexTheta);
    output_tree->Branch("vertexPhi_2cm", &vertexPhi_2cm);
    output_tree->Branch("vertexTheta_2cm", &vertexTheta_2cm);
    output_tree->Branch("vertexPhi_5cm", &vertexPhi_5cm);
    output_tree->Branch("vertexTheta_5cm", &vertexTheta_5cm);
    output_tree->Branch("vertexPhi_7cm", &vertexPhi_7cm);
    output_tree->Branch("vertexTheta_7cm", &vertexTheta_7cm);
    output_tree->Branch("vertexPhi_10cm", &vertexPhi_10cm);
    output_tree->Branch("vertexTheta_10cm", &vertexTheta_10cm);
    output_tree->Branch("vertexPhi_12cm", &vertexPhi_12cm);
    output_tree->Branch("vertexTheta_12cm", &vertexTheta_12cm);
    output_tree->Branch("vertexPhi_15cm", &vertexPhi_15cm);
    output_tree->Branch("vertexTheta_15cm", &vertexTheta_15cm);
    output_tree->Branch("vertexPhi_17cm", &vertexPhi_17cm);
    output_tree->Branch("vertexTheta_17cm", &vertexTheta_17cm);
    output_tree->Branch("vertexPhi_20cm", &vertexPhi_20cm);
    output_tree->Branch("vertexTheta_20cm", &vertexTheta_20cm);
    output_tree->Branch("vertexPhi_30cm", &vertexPhi_30cm);
    output_tree->Branch("vertexTheta_30cm", &vertexTheta_30cm);
    output_tree->Branch("trackAvg15cm_x_sceadded", &trackAvg15cm_x_sceadded);
    output_tree->Branch("trackAvg15cm_y_sceadded", &trackAvg15cm_y_sceadded);
    output_tree->Branch("trackAvg15cm_z_sceadded", &trackAvg15cm_z_sceadded);
    output_tree->Branch("closestWall", &closestWall);
    output_tree->Branch("recoEndPoints_x", &recoEndPoints_x);
    output_tree->Branch("recoEndPoints_y", &recoEndPoints_y);
    output_tree->Branch("recoEndPoints_z", &recoEndPoints_z);
    output_tree->Branch("MCvertex", &MCvertex);
    output_tree->Branch("Ep_t", &Ep_t);
    output_tree->Branch("Em_t", &Em_t);
    output_tree->Branch("Ee_t", &Ee_t);
    output_tree->Branch("MuonStartPoint_X", &MuonStartPoint_X);
    output_tree->Branch("ProtonStartPoint_X", &ProtonStartPoint_X);
    output_tree->Branch("ElectronStartPoint_X", &ElectronStartPoint_X);
    output_tree->Branch("MuonStartPoint_Y", &MuonStartPoint_Y);
    output_tree->Branch("ProtonStartPoint_Y", &ProtonStartPoint_Y);
    output_tree->Branch("ElectronStartPoint_Y", &ElectronStartPoint_Y);
    output_tree->Branch("MuonStartPoint_Z", &MuonStartPoint_Z);
    output_tree->Branch("ProtonStartPoint_Z", &ProtonStartPoint_Z);
    output_tree->Branch("ElectronStartPoint_Z", &ElectronStartPoint_Z);
    output_tree->Branch("MuonEndPoint_X", &MuonEndPoint_X);
    output_tree->Branch("ProtonEndPoint_X", &ProtonEndPoint_X);
    output_tree->Branch("ElectronEndPoint_X", &ElectronEndPoint_X);
    output_tree->Branch("MuonEndPoint_Y", &MuonEndPoint_Y);
    output_tree->Branch("ProtonEndPoint_Y", &ProtonEndPoint_Y);
    output_tree->Branch("ElectronEndPoint_Y", &ElectronEndPoint_Y);
    output_tree->Branch("MuonEndPoint_Z", &MuonEndPoint_Z);
    output_tree->Branch("ProtonEndPoint_Z", &ProtonEndPoint_Z);
    output_tree->Branch("ElectronEndPoint_Z", &ElectronEndPoint_Z);
    
  }
  
}
}
