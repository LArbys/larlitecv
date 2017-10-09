#ifndef __EXTRACT_TRUTH_METHODS_H__
#define __EXTRACT_TRUTH_METHODS_H__

// larlite
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"

class TTree;

namespace larlitecv{

  struct TruthData_t {

    // vertex information
    int mode;
    int current;
    int vertex_boundary;
    float EnuGeV;
    float pos[3];
    float fdwall;

    // secondary information
    int lepton_boundary;       // index of boundary that lepton end point is nearest
    int num_protons_over60mev; // as named
    float dwall_lepton;        // dwll for end of lepton
    float lepton_cosz;         // lepton dir
    float lepton_phiz;         // lepton dir
    float primary_proton_ke;   // ke of the proton driving from the hit nucleon

    void clear();
    void bindToTree( TTree* tree );

  };

  void extractTruth( const larlite::event_mctruth& ev_mctruth, const larlite::event_mctrack& mctrack_v, TruthData_t& data );

  void extractVertexTruth( const larlite::event_mctruth& ev_mctruth, TruthData_t& data );
  void extractLeptonTruth( const larlite::event_mctruth& mctruth_v, const larlite::event_mctrack& mctrack_v, TruthData_t& data );

}

#endif
