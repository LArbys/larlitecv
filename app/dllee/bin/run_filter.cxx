#include "TTree.h"

#include <iostream>

#include "TTree.h"
#include "TFile.h"

#include "DataFormat/IOManager.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/storage_manager.h"

#include "dllee/NumuFilter.h"

int main( int nargs, char** argv ) {

  std::cout << "Run filter" << std::endl;


  // get filter results
  std::map<std::tuple<int,int,int>,bool>     rse_filter;
  std::map<std::tuple<int,int,int,int>,bool> rsev_filter;
  larlitecv::dllee::NumuFilter::ReturnFilteredDictionary( argv[1], rse_filter, rsev_filter );

  std::cout << "[ RSE ] ================" << std::endl;
  for ( auto result=rse_filter.begin(); result!=rse_filter.end(); result++ ) {
    std::cout << "[" << std::get<0>(result->first) << "," << std::get<1>(result->first) << ","
              << std::get<2>(result->first) << "]: "
              << "result=" << result->second << std::endl;
  }
  std::cout << std::endl;
  std::cout << "[ RSE-VERTEXID ] ================" << std::endl;  
  for ( auto result=rsev_filter.begin(); result!=rsev_filter.end(); result++ ) {
    std::cout << "[" << std::get<0>(result->first) << "," << std::get<1>(result->first) << ","
              << std::get<2>(result->first) << "," << std::get<3>(result->first) << "]: "
              << "result=" << result->second << std::endl;
  }
  
  TFile* inputfile  = new TFile( argv[1], "open" );
  TTree* input_tracker_tree  = (TTree*)inputfile->Get("_recoTree_SCEadded");
  TTree* input_vertex_tree   = (TTree*)inputfile->Get("VertexTree");
  TTree* input_eventvtx_tree = (TTree*)inputfile->Get("EventVertexTree");  
  TTree* input_shape_tree    = (TTree*)inputfile->Get("ShapeAnalysis");
  TTree* input_mc_tree       = (TTree*)inputfile->Get("MCTree");
  TTree* input_nufilter_tree = (TTree*)inputfile->Get("NuFilterTree");
  TTree* input_pgraph_tree   = (TTree*)inputfile->Get("PGraphTruthMatch");  
  
  int nentries = input_tracker_tree->GetEntries();
  std::cout << "number of entries: " << nentries << std::endl;

  larcv::IOManager iolcv( larcv::IOManager::kBOTH, "larcv" );
  iolcv.add_in_file( argv[1] );
  iolcv.set_out_file( "test_larcv.root" );
  iolcv.set_verbosity((larcv::msg::Level_t)2);
  iolcv.initialize();
  
  larlite::storage_manager ioll( larlite::storage_manager::kBOTH );
  ioll.add_in_filename( argv[1] );
  ioll.set_out_filename( "test_larlite.root" );
  ioll.set_verbosity((larlite::msg::Level)2);

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

  TTree* output_tracker_tree  = (TTree*)input_tracker_tree->CloneTree(0);
  TTree* output_vertex_tree   = (TTree*)input_vertex_tree->CloneTree(0);
  TTree* output_eventvtx_tree = (TTree*)input_eventvtx_tree->CloneTree(0);
  TTree* output_shape_tree    = (TTree*)input_shape_tree->CloneTree(0);
  TTree* output_mc_tree       = (TTree*)input_mc_tree->CloneTree(0);
  TTree* output_nufilter_tree = (TTree*)input_nufilter_tree->CloneTree(0);
  TTree* output_pgraph_tree   = (TTree*)input_pgraph_tree->CloneTree(0);
  

  // LARCV/LARLITE TREES
  int npass = 0;
  ioll.next_event();  
  for ( int ientry=0; ientry<(int)iolcv.get_n_entries(); ientry++ ) {
    std::cout << "[ENTRY " << ientry << "]" << std::endl;
    iolcv.read_entry(ientry);    
    larcv::EventImage2D* ev_wire = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "wire" );

    std::tuple<int,int,int> rse = { ev_wire->run(), ev_wire->subrun(), ev_wire->event() };
    
    auto it = rse_filter.find( rse );

    if ( it==rse_filter.end() || !it->second )  {
      // skip this entry
      ioll.next_event(false);
      continue;
    }
    
    iolcv.set_id( ev_wire->run(), ev_wire->subrun(), ev_wire->event() );
    iolcv.save_entry();
    ioll.next_event();
    npass++;
  }
  std::cout << "larcv/larlite entries saved: " << npass << std::endl;

  int run, subrun, event, vtxid;
  // VERTEX TREE
  input_vertex_tree->SetBranchAddress("run",    &run);
  input_vertex_tree->SetBranchAddress("subrun", &subrun);
  input_vertex_tree->SetBranchAddress("event",  &event);
  input_vertex_tree->SetBranchAddress("vtxid",  &vtxid);    
  for ( int ientry=0; ientry<input_vertex_tree->GetEntries(); ientry++ ) {
    input_vertex_tree->GetEntry(ientry);
    auto it = rsev_filter.find( std::make_tuple( run, subrun, event, vtxid ) );
    if ( it!=rsev_filter.end() && it->second==true ) {
      output_vertex_tree->Fill();
    }
  }
  
  // TRACKER TREE  
  input_tracker_tree->SetBranchAddress("run",    &run);
  input_tracker_tree->SetBranchAddress("subrun", &subrun);
  input_tracker_tree->SetBranchAddress("event",  &event);
  input_tracker_tree->SetBranchAddress("vtx_id",  &vtxid);    
  for ( int ientry=0; ientry<input_tracker_tree->GetEntries(); ientry++ ) {
    input_tracker_tree->GetEntry(ientry);
    auto it = rsev_filter.find( std::make_tuple( run, subrun, event, vtxid ) );
    if ( it!=rsev_filter.end() && it->second==true ) {
      output_tracker_tree->Fill();
    }
  }

  // SHAPE TREE
  input_shape_tree->SetBranchAddress("run",    &run);
  input_shape_tree->SetBranchAddress("subrun", &subrun);
  input_shape_tree->SetBranchAddress("event",  &event);
  input_shape_tree->SetBranchAddress("vtxid",  &vtxid);    
  for ( int ientry=0; ientry<input_shape_tree->GetEntries(); ientry++ ) {
    input_shape_tree->GetEntry(ientry);
    auto it = rsev_filter.find( std::make_tuple( run, subrun, event, vtxid ) );
    if ( it!=rsev_filter.end() && it->second==true ) {
      output_shape_tree->Fill();
    }
  }

  // PGRAPH TREE  
  input_pgraph_tree->SetBranchAddress("run",    &run);
  input_pgraph_tree->SetBranchAddress("subrun", &subrun);
  input_pgraph_tree->SetBranchAddress("event",  &event);
  input_pgraph_tree->SetBranchAddress("vtxid",  &vtxid);    
  for ( int ientry=0; ientry<input_pgraph_tree->GetEntries(); ientry++ ) {
    input_pgraph_tree->GetEntry(ientry);
    auto it = rsev_filter.find( std::make_tuple( run, subrun, event, vtxid ) );
    if ( it!=rsev_filter.end() && it->second==true ) {
      output_pgraph_tree->Fill();
    }
  }
  
  // EVENT VTX TREE
  input_eventvtx_tree->SetBranchAddress("run",    &run);
  input_eventvtx_tree->SetBranchAddress("subrun", &subrun);
  input_eventvtx_tree->SetBranchAddress("event",  &event);
  for ( int ientry=0; ientry<input_eventvtx_tree->GetEntries(); ientry++ ) {
    input_eventvtx_tree->GetEntry(ientry);
    auto it = rse_filter.find( std::make_tuple( run, subrun, event ) );
    if ( it!=rse_filter.end() && it->second==true ) {
      output_eventvtx_tree->Fill();
    }
  }

  // MC TREE
  input_mc_tree->SetBranchAddress("run",    &run);
  input_mc_tree->SetBranchAddress("subrun", &subrun);
  input_mc_tree->SetBranchAddress("event",  &event);
  for ( int ientry=0; ientry<input_mc_tree->GetEntries(); ientry++ ) {
    input_mc_tree->GetEntry(ientry);
    auto it = rse_filter.find( std::make_tuple( run, subrun, event ) );
    if ( it!=rse_filter.end() && it->second==true ) {
      output_mc_tree->Fill();
    }
  }

  // NUFILTER TREE
  input_nufilter_tree->SetBranchAddress("run",    &run);
  input_nufilter_tree->SetBranchAddress("subrun", &subrun);
  input_nufilter_tree->SetBranchAddress("event",  &event);
  for ( int ientry=0; ientry<input_nufilter_tree->GetEntries(); ientry++ ) {
    input_nufilter_tree->GetEntry(ientry);
    auto it = rse_filter.find( std::make_tuple( run, subrun, event ) );
    if ( it!=rse_filter.end() && it->second==true ) {
      output_nufilter_tree->Fill();
    }
  }

  
  std::cout << "write and close" << std::endl;
  output_eventvtx_tree->Write();
  output_mc_tree->Write();
  output_nufilter_tree->Write();    
  output_tracker_tree->Write();
  output_vertex_tree->Write();
  output_shape_tree->Write();
  output_pgraph_tree->Write();      

  outputfile->Close();
  
  inputfile->Close();
  ioll.close();
  iolcv.finalize();

  return 0;
}
