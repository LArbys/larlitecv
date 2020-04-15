#include "TTree.h"

#include <iostream>

#include "TTree.h"
#include "TFile.h"

#include "DataFormat/IOManager.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/storage_manager.h"
#include "DataFormat/opflash.h"
#include "DataFormat/user_info.h"

#include "dllee/NumuFilter.h"

int main( int nargs, char** argv ) {

  std::cout << "Run filter" << std::endl;

  /*
   * what we are trying to do
   * --------------------------------------------------------------------------------
   * separate trees into event-based and vertex-based
   * we copy the larlite id tree and mc pot tree in full to help with book-keeping
   * read the dlana file and get the finalvertexfile trees
   * make a list of (r,s,e) that passes
   * we pass all events that pass
   * we pass all vertices that pass the event list
   * 
   * we save one file:
   *  one looks like a filtered DLANA file (has both reco and ntuples)
   * --------------------------------------------------------------------------------
   */

  std::string dlmerged_input  = argv[1];
  std::string filtered_output = argv[2];

  /*
  // ALL TREES
  KEY: TTree	image2d_wire_tree;1	wire tree
  KEY: TTree	chstatus_wire_tree;1	wire tree
  KEY: TTree	image2d_ubspurn_plane0_tree;1	ubspurn_plane0 tree
  KEY: TTree	image2d_ubspurn_plane1_tree;1	ubspurn_plane1 tree
  KEY: TTree	image2d_ubspurn_plane2_tree;1	ubspurn_plane2 tree
  KEY: TTree	image2d_thrumu_tree;1	thrumu tree
  KEY: TTree	pgraph_test_tree;1	test tree
  KEY: TTree	pixel2d_test_ctor_tree;1	test_ctor tree
  KEY: TTree	pixel2d_test_img_tree;1	test_img tree
  KEY: TTree	pixel2d_test_super_ctor_tree;1	test_super_ctor tree
  KEY: TTree	pixel2d_test_super_img_tree;1	test_super_img tree
  KEY: TTree	image2d_test_inputimg_tree;1	test_inputimg tree
  KEY: TTree	partroi_croimerge_clip_union_tree;1	croimerge_clip_union tree
  KEY: TTree	image2d_wcshower_tpc_tree;1	wcshower_tpc tree
  KEY: TTree	image2d_wctrack_tpc_tree;1	wctrack_tpc tree
  KEY: TTree	NuFilterTree;1	
  KEY: TTree	MCTree;1	MC infomation
  KEY: TTree	ShapeAnalysis;1	
  KEY: TTree	MatchAnalysis;1	
  KEY: TTree	SecondShowerAnalysis;1	
  KEY: TTree	VertexTree;1	
  KEY: TTree	EventVertexTree;1	
  KEY: TTree	PGraphTruthMatch;1	
  KEY: TTree	_recoTree;1	_recoTree
  KEY: TTree	_recoTree_SCEadded;1	_recoTree_SCEadded
  KEY: TTree	pgraph_inter_par_tree;1	inter_par tree
  KEY: TTree	pixel2d_inter_par_pixel_tree;1	inter_par_pixel tree
  KEY: TTree	pixel2d_inter_img_pixel_tree;1	inter_img_pixel tree
  KEY: TTree	pixel2d_inter_int_pixel_tree;1	inter_int_pixel tree
  KEY: TTree	SelNueID;1	
  KEY: TTree	larlite_id_tree;1	LArLite Event ID Tree
  KEY: TTree	daqheadertimeuboone_daq_tree;1	daqheadertimeuboone Tree by daq
  KEY: TTree	hit_gaushit_tree;1	hit Tree by gaushit
  KEY: TTree	hit_portedThresholdhit_tree;1	hit Tree by portedThresholdhit
  KEY: TTree	crthit_crthitcorr_tree;1	crthit Tree by crthitcorr
  KEY: TTree	crttrack_crttrack_tree;1	crttrack Tree by crttrack
  KEY: TTree	ophit_ophitBeam_tree;1	ophit Tree by ophitBeam
  KEY: TTree	ophit_ophitBeamCalib_tree;1	ophit Tree by ophitBeamCalib
  KEY: TTree	ophit_ophitCosmic_tree;1	ophit Tree by ophitCosmic
  KEY: TTree	ophit_ophitCosmicCalib_tree;1	ophit Tree by ophitCosmicCalib
  KEY: TTree	opflash_opflashBeam_tree;1	opflash Tree by opflashBeam
  KEY: TTree	opflash_opflashCosmic_tree;1	opflash Tree by opflashCosmic
  KEY: TTree	opflash_portedFlash_tree;1	opflash Tree by portedFlash
  KEY: TTree	opflash_simpleFlashBeam_tree;1	opflash Tree by simpleFlashBeam
  KEY: TTree	opflash_simpleFlashBeam::DLWCDeploy_tree;1	opflash Tree by simpleFlashBeam::DLWCDeploy
  KEY: TTree	opflash_simpleFlashCosmic_tree;1	opflash Tree by simpleFlashCosmic
  KEY: TTree	opflash_simpleFlashCosmic::DLWCDeploy_tree;1	opflash Tree by simpleFlashCosmic::DLWCDeploy
  KEY: TTree	sps_portedSpacePointsThreshold_tree;1	sps Tree by portedSpacePointsThreshold
  KEY: TTree	track_inter_track_tree;1	track Tree by inter_track
  KEY: TTree	track_trackReco_tree;1	track Tree by trackReco
  KEY: TTree	track_trackReco_sceadded_tree;1	track Tree by trackReco_sceadded
  KEY: TTree	shower_ssnetshowerreco_tree;1	shower Tree by ssnetshowerreco
  KEY: TTree	vertex_inter_vertex_tree;1	vertex Tree by inter_vertex
  KEY: TTree	vertex_trackReco_tree;1	vertex Tree by trackReco
  KEY: TTree	trigger_daq_tree;1	trigger Tree by daq
  KEY: TTree	ass_inter_ass_tree;1	ass Tree by inter_ass
  KEY: TTree	ass_opflashBeam_tree;1	ass Tree by opflashBeam
  KEY: TTree	ass_opflashCosmic_tree;1	ass Tree by opflashCosmic
  KEY: TTree	ass_portedFlash_tree;1	ass Tree by portedFlash
  KEY: TTree	ass_portedSpacePointsThreshold_tree;1	ass Tree by portedSpacePointsThreshold
  KEY: TTree	ass_simpleFlashBeam_tree;1	ass Tree by simpleFlashBeam
  KEY: TTree	ass_simpleFlashBeam::DLWCDeploy_tree;1	ass Tree by simpleFlashBeam::DLWCDeploy
  KEY: TTree	ass_simpleFlashCosmic_tree;1	ass Tree by simpleFlashCosmic
  KEY: TTree	ass_simpleFlashCosmic::DLWCDeploy_tree;1	ass Tree by simpleFlashCosmic::DLWCDeploy
  KEY: TTree	ass_trackReco_tree;1	ass Tree by trackReco
  KEY: TTree	ass_trackReco_sceadded_tree;1	ass Tree by trackReco_sceadded
  KEY: TTree	swtrigger_swtrigger_tree;1	swtrigger Tree by swtrigger
  KEY: TTree	larflowcluster_ssnetshowerreco_tree;1	larflowcluster Tree by ssnetshowerreco
  KEY: TTree	clustermask_mrcnn_masks_tree;1	mrcnn_masks tree
  KEY: TTree	sparseimg_larflow_tree;1	larflow tree
  KEY: TTree	sparseimg_sparseuresnetout_tree;1	sparseuresnetout tree
  KEY: TTree	shower_ssnetshowerrecov2ana_tree;1	shower Tree by ssnetshowerrecov2ana
  KEY: TTree	shower_ssnetshowerrecov2ana_sec_tree;1	shower Tree by ssnetshowerrecov2ana_sec
  KEY: TTree	larflowcluster_ssnetshowerrecov2ana_tree;1	larflowcluster Tree by ssnetshowerrecov2ana
  KEY: TTree	track_dqdx_U_tree;1	track Tree by dqdx_U
  KEY: TTree	track_dqdx_V_tree;1	track Tree by dqdx_V
  KEY: TTree	track_dqdx_Y_tree;1	track Tree by dqdx_Y
  KEY: TDirectoryFile	dlana;1	dlana
  KEY: TDirectoryFile	mpid;1	mpid
  KEY: TDirectoryFile	ssnetshowerreco;1	ssnetshowerreco
  KEY: TH2D	residual_dqdx;1	residual_dqdx 
  */

  
  // list of event-based trees
  std::vector<std::string> event_indexed_trees = 
    {
      "larlite_id_tree",
      "image2d_wire_tree",
      "chstatus_wire_tree",
      "image2d_ubspurn_plane0_tree",
      "image2d_ubspurn_plane1_tree",
      "image2d_ubspurn_plane2_tree",
      "image2d_thrumu_tree",
      "pgraph_test_tree",
      "pixel2d_test_ctor_tree",
      "pixel2d_test_img_tree",
      "pixel2d_test_super_ctor_tree",
      "pixel2d_test_super_img_tree",
      "image2d_test_inputimg_tree",
      "partroi_croimerge_clip_union_tree",
      //"image2d_wcshower_tpc_tree",
      //"image2d_wctrack_tpc_tree",
      "NuFilterTree",
      "MCTree",
      "EventVertexTree",
      "PGraphTruthMatch",      
      "pgraph_inter_par_tree",
      "pixel2d_inter_par_pixel_tree",
      "pixel2d_inter_img_pixel_tree",
      "pixel2d_inter_int_pixel_tree",
      "daqheadertimeuboone_daq_tree",
      "hit_gaushit_tree",
      "hit_portedThresholdhit_tree",
      "crthit_crthitcorr_tree",
      "crttrack_crttrack_tree",
      "ophit_ophitBeam_tree",
      "ophit_ophitBeamCalib_tree",
      "ophit_ophitCosmic_tree",
      "ophit_ophitCosmicCalib_tree",
      "opflash_opflashBeam_tree",
      "opflash_opflashCosmic_tree",
      "opflash_portedFlash_tree",
      "opflash_simpleFlashBeam_tree",
      "opflash_simpleFlashBeam::DLWCDeploy_tree",
      "opflash_simpleFlashCosmic_tree",
      "opflash_simpleFlashCosmic::DLWCDeploy_tree",
      "sps_portedSpacePointsThreshold_tree",
      "track_inter_track_tree",
      "track_trackReco_tree",
      "track_trackReco_sceadded_tree",
      "shower_ssnetshowerreco_tree",
      "vertex_inter_vertex_tree",
      "vertex_trackReco_tree",
      "trigger_daq_tree",
      "ass_inter_ass_tree",
      "ass_opflashBeam_tree",
      "ass_opflashCosmic_tree",
      "ass_portedFlash_tree",
      "ass_portedSpacePointsThreshold_tree",
      "ass_simpleFlashBeam_tree",
      "ass_simpleFlashBeam::DLWCDeploy_tree",
      "ass_simpleFlashCosmic_tree",
      "ass_simpleFlashCosmic::DLWCDeploy_tree",
      "ass_trackReco_tree",
      "ass_trackReco_sceadded_tree",
      "swtrigger_swtrigger_tree",
      "larflowcluster_ssnetshowerreco_tree",
      "clustermask_mrcnn_masks_tree",
      "sparseimg_larflow_tree",
      "sparseimg_sparseuresnetout_tree",
      "shower_ssnetshowerrecov2ana_tree",
      "shower_ssnetshowerrecov2ana_sec_tree",
      "larflowcluster_ssnetshowerrecov2ana_tree",
      "track_dqdx_U_tree",
      "track_dqdx_V_tree",
      "track_dqdx_Y_tree",
      "ssnetshowerreco/ssnetshowerrecov2ana_anatree"
    };
  

  // list of vertex-based trees
  std::vector<std::string> vertex_indexed_trees =
    {
      "VertexTree",
      "ShapeAnalysis",
      "MatchAnalysis",
      "SecondShowerAnalysis",
      "_recoTree",
      "_recoTree_SCEadded",
      "SelNueID",
      "mpid/multipid_tree",
      "dlana/FinalVertexVariables"
    };
  // list of copy-in-full trees

  TFile* inputfile  = new TFile( dlmerged_input.c_str(), "open" );

  int num_event_entries = 0;
  std::vector< TTree* > event_indexed_trees_v;
  for ( auto const& name : event_indexed_trees ) {
    event_indexed_trees_v.push_back( (TTree*)inputfile->Get(name.c_str()) );

    if ( event_indexed_trees_v.back()==NULL  ) {
      std::cout << "Error loading event-indexed tree: " << name << std::endl;
    }
    
    if ( event_indexed_trees_v.size()==1 ) {
      num_event_entries = event_indexed_trees_v.back()->GetEntries();
      std::cout << "Event-indexed trees, num entries=" << num_event_entries << std::endl;
    }
    else {
      if ( num_event_entries!=event_indexed_trees_v.back()->GetEntries() ) {
        std::cout << "Tree [" << name << "] does not match the entry-indexed tree entries"  << std::endl;
      }
    }
  }
  std::cout << "Number of event-indexed trees: " << event_indexed_trees_v.size() << std::endl;

  int num_vertex_entries = 0;
  std::vector< TTree* > vertex_indexed_trees_v;
  for ( auto const& name : vertex_indexed_trees ) {
    vertex_indexed_trees_v.push_back( (TTree*)inputfile->Get(name.c_str()) );

    if ( vertex_indexed_trees_v.back()==NULL  ) {
      std::cout << "Error loading vertex-indexed tree: " << name << std::endl;
    }
    
    if ( vertex_indexed_trees_v.size()==1 ) {
      num_vertex_entries = vertex_indexed_trees_v.back()->GetEntries();
      std::cout << "Vertex-indexed trees, num entries=" << num_vertex_entries << std::endl;
    }
    else {
      if ( num_vertex_entries!=vertex_indexed_trees_v.back()->GetEntries() ) {
        std::cout << "Tree [" << name << "] does not match the entry-indexed tree entries"  << std::endl;
      }
    }
  }
  std::cout << "Number of vertex-indexed trees: " << vertex_indexed_trees_v.size() << std::endl;  

  // get filter results
  std::map<std::tuple<int,int,int>,bool>     rse_filter;
  std::map<std::tuple<int,int,int,int>,bool> rsev_filter;
  larlitecv::dllee::NumuFilter::ReturnFilteredDictionary( dlmerged_input, rse_filter, rsev_filter, true );

  // save tree entries
  TFile* outputfile = new TFile( filtered_output.c_str(), "new");

  // save entry-indexed trees
  std::vector<TTree*> out_event_indexed_v;
  for ( auto& ptree : event_indexed_trees_v ) {
    TTree* outtree = ptree->CloneTree(0);
    out_event_indexed_v.push_back(outtree);
  }
  
  unsigned int run;
  unsigned int subrun;
  unsigned int event;
  event_indexed_trees_v[0]->SetBranchAddress( "_run_id",    &run );
  event_indexed_trees_v[0]->SetBranchAddress( "_subrun_id", &subrun );
  event_indexed_trees_v[0]->SetBranchAddress( "_event_id",  &event );  

  for (int ientry=0; ientry<num_event_entries; ientry++ ) {

    event_indexed_trees_v[0]->GetEntry(ientry);
    
    std::tuple<int,int,int> rse = std::make_tuple(run,subrun,event);
    auto it = rse_filter.find(rse);
    if ( it!=rse_filter.end() && it->second ) {
      // passed. save event entry.    
      for ( int itree = 0; itree<(int)event_indexed_trees_v.size(); itree++ ) {
        event_indexed_trees_v[itree]->GetEntry(ientry);
        out_event_indexed_v[itree]->Fill();
      }
    }
  }

  // save vertex-indexed trees
  int vtx_run;
  int vtx_subrun;
  int vtx_event;
  int vtx_vtxid;
  vertex_indexed_trees_v[0]->SetBranchAddress("run",    &vtx_run);
  vertex_indexed_trees_v[0]->SetBranchAddress("subrun", &vtx_subrun);
  vertex_indexed_trees_v[0]->SetBranchAddress("event",  &vtx_event);
  vertex_indexed_trees_v[0]->SetBranchAddress("vtxid",  &vtx_vtxid);
  
  std::vector<TTree*> out_vertex_indexed_v;
  for ( auto& ptree : vertex_indexed_trees_v ) {
    TTree* outtree = ptree->CloneTree(0);
    out_vertex_indexed_v.push_back(outtree);
  }
  
  for (int ientry=0; ientry<num_vertex_entries; ientry++ ) {

    vertex_indexed_trees_v[0]->GetEntry(ientry);
    
    std::tuple<int,int,int,int> rse = std::make_tuple(vtx_run,vtx_subrun,vtx_event,vtx_vtxid);
    auto it = rsev_filter.find(rse);
    if ( it!=rsev_filter.end() && it->second ) {
      // passed. save event entry.    
      for ( int itree = 0; itree<(int)vertex_indexed_trees_v.size(); itree++ ) {
        vertex_indexed_trees_v[itree]->GetEntry(ientry);
        out_vertex_indexed_v[itree]->Fill();
      }
    }
  }
  
  for ( auto& ptree : out_event_indexed_v )
    ptree->Write();
  for ( auto& ptree : out_vertex_indexed_v )
    ptree->Write();

  outputfile->Close();
  inputfile->Close();
  //ioll.close();
  //iolcv.finalize();

  return 0;
}
