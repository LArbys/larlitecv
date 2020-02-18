#include "TTree.h"

#include <iostream>

#include "TTree.h"
#include "TFile.h"

#include "DataFormat/IOManager.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/storage_manager.h"

#include "dllee/TrackRecoFinalFileVariables.h"


int main( int nargs, char** argv ) {

  std::cout << "Run filter" << std::endl;

  TFile* inputfile  = new TFile( argv[1], "open" );
  TTree* input_tree  = (TTree*)inputfile->Get("_recoTree_SCEadded");

  int nentries = input_tree->GetEntries();
  std::cout << "number of entries: " << nentries << std::endl;

  larcv::IOManager iolcv( larcv::IOManager::kBOTH, "larcv" );
  iolcv.add_in_file( argv[1] );
  iolcv.set_out_file( "test_larcv.root" );
  iolcv.set_verbosity((larcv::msg::Level_t)1);
  iolcv.initialize();
  
  larlite::storage_manager ioll( larlite::storage_manager::kBOTH );
  ioll.add_in_filename( argv[1] );
  ioll.set_out_filename( "test_larlite.root" );
  ioll.set_verbosity((larlite::msg::Level)1);

  ioll.set_data_to_read( larlite::data::kDAQHeaderTimeUBooNE, "daq" );
  ioll.set_data_to_read( larlite::data::kDAQHeaderTimeUBooNE, "triggersim" );  
  ioll.set_data_to_read( larlite::data::kHit, "dl" );
  ioll.set_data_to_read( larlite::data::kHit, "dlrea" );
  ioll.set_data_to_read( larlite::data::kHit, "gaushit" );  
  ioll.set_data_to_read( larlite::data::kCRTHit, "crthitcorr" );
  ioll.set_data_to_read( larlite::data::kCRTTrack, "crttrack" );
  ioll.set_data_to_read( larlite::data::kOpHit, "ophitBeam" );
  ioll.set_data_to_read( larlite::data::kOpHit, "ophitCosmic" );
  ioll.set_data_to_read( larlite::data::kOpFlash, "opflash" );
  ioll.set_data_to_read( larlite::data::kOpFlash, "opflashBeam" );
  ioll.set_data_to_read( larlite::data::kOpFlash, "opflashCosmic" );      
  ioll.set_data_to_read( larlite::data::kOpFlash, "ophypo" );
  ioll.set_data_to_read( larlite::data::kOpFlash, "simpleFlashBeam" );
  ioll.set_data_to_read( larlite::data::kOpFlash, "simpleFlashCosmic" );
  ioll.set_data_to_read( larlite::data::kCluster, "dl" );
  ioll.set_data_to_read( larlite::data::kCluster, "dlraw" );
  ioll.set_data_to_read( larlite::data::kTrack,   "all3d" );
  ioll.set_data_to_read( larlite::data::kTrack,   "croi3d" );
  ioll.set_data_to_read( larlite::data::kTrack,   "dl" );
  ioll.set_data_to_read( larlite::data::kTrack,   "inter_track" );
  ioll.set_data_to_read( larlite::data::kTrack,   "mergedstopmu3d" );
  ioll.set_data_to_read( larlite::data::kTrack,   "mergedthrumu3d" );
  ioll.set_data_to_read( larlite::data::kTrack,   "mergeduntagged3d" );
  ioll.set_data_to_read( larlite::data::kTrack,   "stopmu3d" );
  ioll.set_data_to_read( larlite::data::kTrack,   "streclustered3d" );
  ioll.set_data_to_read( larlite::data::kTrack,   "thrumu3d" );
  ioll.set_data_to_read( larlite::data::kTrack,   "trackReco" );
  ioll.set_data_to_read( larlite::data::kTrack,   "trackReco_sceadded" );
  ioll.set_data_to_read( larlite::data::kTrack,   "untagged3d" );
  ioll.set_data_to_read( larlite::data::kShower,  "dl" );    
  ioll.set_data_to_read( larlite::data::kShower,  "showerreco" );
  ioll.set_data_to_read( larlite::data::kShower,  "ssnetshowerreco" );              
  ioll.set_data_to_read( larlite::data::kVertex,  "dl" );
  ioll.set_data_to_read( larlite::data::kVertex,  "dlraw" );
  ioll.set_data_to_read( larlite::data::kVertex,  "inter_vertex" );
  ioll.set_data_to_read( larlite::data::kVertex,  "trackReco" );
  ioll.set_data_to_read( larlite::data::kPFParticle,  "dl" );
  ioll.set_data_to_read( larlite::data::kPFParticle,  "dlraw" );
  ioll.set_data_to_read( larlite::data::kUserInfo,    "croicutresults" );
  ioll.set_data_to_read( larlite::data::kUserInfo,    "endpointresults" );
  ioll.set_data_to_read( larlite::data::kUserInfo,    "precutresults" );
  ioll.set_data_to_read( larlite::data::kTrigger,     "daq" );
  ioll.set_data_to_read( larlite::data::kTrigger,     "triggersim" );        
  ioll.set_data_to_read( larlite::data::kAssociation, "dl" );
  ioll.set_data_to_read( larlite::data::kAssociation, "dlraw" );
  ioll.set_data_to_read( larlite::data::kAssociation, "inter_ass" );
  ioll.set_data_to_read( larlite::data::kAssociation, "opflashBeam" );
  ioll.set_data_to_read( larlite::data::kAssociation, "opflashCosmic" );
  ioll.set_data_to_read( larlite::data::kAssociation, "showerreco" );
  ioll.set_data_to_read( larlite::data::kAssociation, "simpleFlashBeam" );
  ioll.set_data_to_read( larlite::data::kAssociation, "simpleFlashCosmic" );
  ioll.set_data_to_read( larlite::data::kAssociation, "trackReco" );
  ioll.set_data_to_read( larlite::data::kAssociation, "trackReco_sceaddded" );
  ioll.set_data_to_read( larlite::data::kSWTrigger,   "swtrigger" );      
  ioll.open();

  TFile* outputfile = new TFile( argv[2], "new");  
  TTree* output_tree = new TTree("_recoTree_SCEadded","filterd Track Reco tree");

  larlitecv::dllee::TrackRecoFinalFileVariables( input_tree, output_tree );  


  ioll.next_event();  
  for ( int ientry=0; ientry<(int)iolcv.get_n_entries(); ientry++ ) {
    std::cout << "[ENTRY " << ientry << "]" << std::endl;
    iolcv.read_entry(ientry);    
    larcv::EventImage2D* ev_wire = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "wire" );
    
    input_tree->GetEntry(ientry);

    output_tree->Fill();
    
    iolcv.set_id( ev_wire->run(), ev_wire->subrun(), ev_wire->event() );
    iolcv.save_entry();
    ioll.next_event();    
  }

  std::cout << "write and close" << std::endl;
  output_tree->Write();
  outputfile->Close();
  
  inputfile->Close();
  ioll.close();
  iolcv.finalize();

  return 0;
}
