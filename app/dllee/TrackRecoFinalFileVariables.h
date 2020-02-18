#ifndef __DLLEE_TRACKRECO_FINALFILE_VARIABLES_H__
#define __DLLEE_TRACKRECO_FINALFILE_VARIABLES_H__

#include <vector>
#include "TTree.h"
#include "TBranch.h"
#include "TVector3.h"

using namespace std;

namespace larlitecv {
namespace dllee {

  class TrackRecoFinalFileVariables {

  public:

    TrackRecoFinalFileVariables( TTree* input_tree, TTree* output_tree );
    virtual ~TrackRecoFinalFileVariables();

  protected:
    
    virtual void defineInputTreeVariables( TTree* );
    virtual void defineOutputTreeVariables( TTree* );

  public:

    void clearVariables();

    // Declaration of leaf types
    Int_t           run;
    Int_t           subrun;
    Int_t           event;
    Int_t           entry;
    Int_t           vtx_id;
    Float_t         vtx_x;
    Float_t         vtx_y;
    Float_t         vtx_z;
    Int_t           randomSeed;
    vector<int>     *trk_id_v;
    Int_t           NtracksReco;
    Double_t        MinLengthAllowed;
    vector<vector<double> > *trackQ3_v;
    vector<vector<double> > *trackQ5_v;
    vector<vector<double> > *trackQ10_v;
    vector<vector<double> > *trackQ20_v;
    vector<vector<double> > *trackQ30_v;
    vector<vector<double> > *trackQ50_v;
    vector<double>  *E_muon_v;
    vector<double>  *E_proton_v;
    vector<double>  *Length_v;
    vector<double>  *Avg_Ion_v;
    vector<double>  *Avg_IonY_v;
    vector<double>  *Ion_5cm_v;
    vector<double>  *Ion_10cm_v;
    vector<double>  *Ion_tot_v;
    vector<double>  *IonY_5cm_v;
    vector<double>  *IonY_10cm_v;
    vector<double>  *IonY_tot_v;
    vector<double>  *Truncated_dQdX1_v;
    vector<double>  *Truncated_dQdX3_v;
    vector<double>  *IondivLength_v;
    vector<vector<double> > *TotalADCvalues_v;
    vector<vector<double> > *Angle_v;
    vector<int>     *Reco_goodness_v;
    vector<int>     *track_Goodness_v;
    Bool_t          GoodVertex;
    Bool_t          missingTrack;
    Bool_t          nothingReconstructed;
    Bool_t          tooShortDeadWire;
    Bool_t          tooShortFaintTrack;
    Bool_t          tooManyTracksAtVertex;
    Bool_t          possibleCosmic;
    Bool_t          possiblyCrossing;
    Bool_t          branchingTracks;
    Bool_t          jumpingTracks;
    TVector3        *RecoVertex;
    vector<double>  *vertexPhi;
    vector<double>  *vertexTheta;
    vector<double>  *vertexPhi_2cm;
    vector<double>  *vertexTheta_2cm;
    vector<double>  *vertexPhi_5cm;
    vector<double>  *vertexTheta_5cm;
    vector<double>  *vertexPhi_7cm;
    vector<double>  *vertexTheta_7cm;
    vector<double>  *vertexPhi_10cm;
    vector<double>  *vertexTheta_10cm;
    vector<double>  *vertexPhi_12cm;
    vector<double>  *vertexTheta_12cm;
    vector<double>  *vertexPhi_15cm;
    vector<double>  *vertexTheta_15cm;
    vector<double>  *vertexPhi_17cm;
    vector<double>  *vertexTheta_17cm;
    vector<double>  *vertexPhi_20cm;
    vector<double>  *vertexTheta_20cm;
    vector<double>  *vertexPhi_30cm;
    vector<double>  *vertexTheta_30cm;
    vector<double>  *trackAvg15cm_x_sceadded;
    vector<double>  *trackAvg15cm_y_sceadded;
    vector<double>  *trackAvg15cm_z_sceadded;
    vector<double>  *closestWall;
    vector<double>  *recoEndPoints_x;
    vector<double>  *recoEndPoints_y;
    vector<double>  *recoEndPoints_z;
    TVector3        *MCvertex;
    Double_t        Ep_t;
    Double_t        Em_t;
    Double_t        Ee_t;
    Double_t        MuonStartPoint_X;
    Double_t        ProtonStartPoint_X;
    Double_t        ElectronStartPoint_X;
    Double_t        MuonStartPoint_Y;
    Double_t        ProtonStartPoint_Y;
    Double_t        ElectronStartPoint_Y;
    Double_t        MuonStartPoint_Z;
    Double_t        ProtonStartPoint_Z;
    Double_t        ElectronStartPoint_Z;
    Double_t        MuonEndPoint_X;
    Double_t        ProtonEndPoint_X;
    Double_t        ElectronEndPoint_X;
    Double_t        MuonEndPoint_Y;
    Double_t        ProtonEndPoint_Y;
    Double_t        ElectronEndPoint_Y;
    Double_t        MuonEndPoint_Z;
    Double_t        ProtonEndPoint_Z;
    Double_t        ElectronEndPoint_Z;

    // List of branches
    TBranch        *b__run;   //!
    TBranch        *b__subrun;   //!
    TBranch        *b__event;   //!
    TBranch        *b__entry;   //!
    TBranch        *b_vtx_id;   //!
    TBranch        *b_vtx_x_sceadded;   //!
    TBranch        *b_vtx_y_sceadded;   //!
    TBranch        *b_vtx_z_sceadded;   //!
    TBranch        *b_randomSeed;   //!
    TBranch        *b_trk_id_v;   //!
    TBranch        *b_NtracksReco;   //!
    TBranch        *b_MinLengthAllowed;   //!
    TBranch        *b_trackQ3_v;   //!
    TBranch        *b_trackQ5_v;   //!
    TBranch        *b_trackQ10_v;   //!
    TBranch        *b_trackQ20_v;   //!
    TBranch        *b_trackQ30_v;   //!
    TBranch        *b_trackQ50_v;   //!
    TBranch        *b_E_muon_v;   //!
    TBranch        *b_E_proton_v;   //!
    TBranch        *b_Length_v;   //!
    TBranch        *b_Avg_Ion_v;   //!
    TBranch        *b_Avg_IonY_v;   //!
    TBranch        *b_Ion_5cm_v;   //!
    TBranch        *b_Ion_10cm_v;   //!
    TBranch        *b_Ion_tot_v;   //!
    TBranch        *b_IonY_5cm_v;   //!
    TBranch        *b_IonY_10cm_v;   //!
    TBranch        *b_IonY_tot_v;   //!
    TBranch        *b_Truncated_dQdX1_v;   //!
    TBranch        *b_Truncated_dQdX3_v;   //!
    TBranch        *b_IondivLength_v;   //!
    TBranch        *b_TotalADCvalues_v;   //!
    TBranch        *b_Angle_v;   //!
    TBranch        *b_Reco_goodness_v;   //!
    TBranch        *b_track_Goodness_v;   //!
    TBranch        *b_GoodVertex;   //!
    TBranch        *b_missingTrack;   //!
    TBranch        *b_nothingReconstructed;   //!
    TBranch        *b_tooShortDeadWire;   //!
    TBranch        *b_tooShortFaintTrack;   //!
    TBranch        *b_tooManyTracksAtVertex;   //!
    TBranch        *b_possibleCosmic;   //!
    TBranch        *b_possiblyCrossing;   //!
    TBranch        *b_branchingTracks;   //!
    TBranch        *b_jumpingTracks;   //!
    TBranch        *b_RecoVertex;   //!
    TBranch        *b_vertexPhi;   //!
    TBranch        *b_vertexTheta;   //!
    TBranch        *b_vertexPhi_2cm;   //!
    TBranch        *b_vertexTheta_2cm;   //!
    TBranch        *b_vertexPhi_5cm;   //!
    TBranch        *b_vertexTheta_5cm;   //!
    TBranch        *b_vertexPhi_7cm;   //!
    TBranch        *b_vertexTheta_7cm;   //!
    TBranch        *b_vertexPhi_10cm;   //!
    TBranch        *b_vertexTheta_10cm;   //!
    TBranch        *b_vertexPhi_12cm;   //!
    TBranch        *b_vertexTheta_12cm;   //!
    TBranch        *b_vertexPhi_15cm;   //!
    TBranch        *b_vertexTheta_15cm;   //!
    TBranch        *b_vertexPhi_17cm;   //!
    TBranch        *b_vertexTheta_17cm;   //!
    TBranch        *b_vertexPhi_20cm;   //!
    TBranch        *b_vertexTheta_20cm;   //!
    TBranch        *b_vertexPhi_30cm;   //!
    TBranch        *b_vertexTheta_30cm;   //!
    TBranch        *b_trackAvg15cm_x_sceadded;   //!
    TBranch        *b_trackAvg15cm_y_sceadded;   //!
    TBranch        *b_trackAvg15cm_z_sceadded;   //!
    TBranch        *b_closestWall;   //!
    TBranch        *b_recoEndPoints_x;   //!
    TBranch        *b_recoEndPoints_y;   //!
    TBranch        *b_recoEndPoints_z;   //!
    TBranch        *b_MCvertex;   //!
    TBranch        *b_Ep_t;   //!
    TBranch        *b_Em_t;   //!
    TBranch        *b_Ee_t;   //!
    TBranch        *b_MuonStartPoint_X;   //!
    TBranch        *b_ProtonStartPoint_X;   //!
    TBranch        *b_ElectronStartPoint_X;   //!
    TBranch        *b_MuonStartPoint_Y;   //!
    TBranch        *b_ProtonStartPoint_Y;   //!
    TBranch        *b_ElectronStartPoint_Y;   //!
    TBranch        *b_MuonStartPoint_Z;   //!
    TBranch        *b_ProtonStartPoint_Z;   //!
    TBranch        *b_ElectronStartPoint_Z;   //!
    TBranch        *b_MuonEndPoint_X;   //!
    TBranch        *b_ProtonEndPoint_X;   //!
    TBranch        *b_ElectronEndPoint_X;   //!
    TBranch        *b_MuonEndPoint_Y;   //!
    TBranch        *b_ProtonEndPoint_Y;   //!
    TBranch        *b_ElectronEndPoint_Y;   //!
    TBranch        *b_MuonEndPoint_Z;   //!
    TBranch        *b_ProtonEndPoint_Z;   //!
    TBranch        *b_ElectronEndPoint_Z;   //!

  };
  
}
}

#endif
