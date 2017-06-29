#ifndef __CosmicTagger_h__
#define __CosmicTagger_h__

/* ================================================
 * CosmicTagger
 * 
 * This class encapulates the tagger code
 * ================================================= */

#include <vector>
#include <string>

// larlite

// larcv
#include "Base/PSet.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "TaggerCROIAlgoConfig.h"
#include "TaggerCROITypes.h"

namespace larlitecv {

  class EmptyChannelAlgo;
  class UnipolarHackAlgo;
  class TaggerCROIAlgo;
  
  class CosmicTagger {

  protected:
    CosmicTagger();
    
  public:
    CosmicTagger( std::string tagger_cfg );
    CosmicTagger( std::string tagger_cfg, std::string larcv_list, std::string larlite_list );
    CosmicTagger( std::string tagger_cfg, const std::vector<std::string>& larcv_filepaths, const std::vector<std::string>& larlite_filepaths );

    virtual ~CosmicTagger();

    void setEntry( int ientry ) { m_entry = ientry; };
    void runOneEvent( int ientry );
    virtual bool processInputImages();
    bool findBoundaryEnds();
    bool findThruMuTracks();
    bool runCombinedThruMu(); //< run boundary endpoint finder and thrumu tracker in one step    
    bool findStopMu();
    //bool findUntaggedClusters();
    bool findCROI();
	  	
    void PrintTruthVertexInfo();	 
	
    larlitecv::TaggerCROIAlgoConfig& getConfig() { return m_tagger_cfg; }; //< can change any parameter here

    std::string printState();

    bool isInputReady() { return m_state.input_ready; };
    bool hasBoundaryTaggerRun() { return m_state.boundary_run; };    
    bool hasThruMuTrackerRun() { return m_state.thrumu_run; };
    bool hasStopMuRun() { return m_state.stopmu_run; };
    //bool hasUntaggedRun() { return m_state.untagged_run; };
    bool hasCROIrun() { return m_state.croi_run; };

    const larlitecv::InputPayload&  getInputData()  { return m_input_data; };
    const larlitecv::ThruMuPayload& getThruMuData() { return m_thrumu_data; };
    const larlitecv::StopMuPayload& getStopMuData() { return m_stopmu_data; };
    const larlitecv::CROIPayload&   getCROIData()   { return m_croi_data_v.front(); };
    
  protected:

    larlitecv::DataCoordinator m_dataco_input;
    larlitecv::DataCoordinator m_dataco_output;
    int m_entry;
    int m_nentries;
    
    std::string m_larcv_image_producer;
    std::string m_larcv_chstatus_producer;
    bool m_DeJebWires;
    float m_jebwiresfactor;
    std::vector<float> m_emptych_thresh;
    std::string m_chstatus_datatype;
    std::vector<std::string> m_opflash_producers;
    bool m_RunThruMu;
    bool m_RunStopMu;
    bool m_RunCROI;
    bool m_save_thrumu_space;
    bool m_save_stopmu_space;
    bool m_save_croi_space;
    bool m_save_mc;
    bool m_skip_empty_events;
    bool m_apply_unipolar_hack;
    void setRunParameters();

    // Config file
    larcv::PSet m_toplevel_pset;
    larcv::PSet m_pset;
    std::string m_cfg_file;
    larlitecv::TaggerCROIAlgoConfig m_tagger_cfg;
    void setConfigFile( std::string cfgfile );

    larlitecv::EmptyChannelAlgo* m_emptyalgo;
    larlitecv::UnipolarHackAlgo* m_unihackalgo;
    larlitecv::TaggerCROIAlgo*   m_taggercroialgo;
    void configure_algos();    

    // We track what has run by this state. In principle we can clear downstream states/results and rerun.
    struct State_t {
      bool configured;
      bool input_ready;
      bool boundary_run;
      bool thrumu_run;
      bool stopmu_run;
      bool untagged_run;
      bool croi_run;
      State_t() {
	clear();
      };
      void clear() {
	configured = input_ready = boundary_run = thrumu_run = stopmu_run = untagged_run = croi_run = false;
      };
    } m_state;

    // Data Classes
    larlitecv::InputPayload  m_input_data;
    larlitecv::ThruMuPayload m_thrumu_data;
    larlitecv::StopMuPayload m_stopmu_data;
    std::vector<larlitecv::CROIPayload> m_croi_data_v;
   
  };
  
}

#endif
