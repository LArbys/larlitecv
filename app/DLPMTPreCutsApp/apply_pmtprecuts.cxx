#include <iostream>
#include <string>

#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

#include "larlite/UserDev/BasicTool/FhiclLite/PSet.h"
#include "larlite/UserDev/SelectionTool/LEEPreCuts/LEEPreCut.h"

int main( int nargs, char** argv ) {

  // pre-amble/parse arguments
  std::cout << "Apply PMT Pre-Cuts" << std::endl;

  if ( nargs!=2 ) {
    std::cout << "usage: ./apply_pmtprecuts [config file]" << std::endl;
    return 0;
  }

  std::string cfg_file = argv[1];

  // tagger routine configuration
  larcv::PSet cfg  = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg.get<larcv::PSet>("ApplyPMTPrecuts");
  std::string larcv_inputlist   = pset.get<std::string>("LArCVInputFilelist",   "larcv" );
  std::string larlite_inputlist = pset.get<std::string>("LArLiteInputFilelist", "larlite" );  
  std::string ophit_producer    = pset.get<std::string>("OpHitProducer");

  // Setup PMT Precuts algo
  fcllite::PSet tmp( "tmp",pset.get<larcv::PSet>("LEEPreCut").data_string() ); 
  fcllite::PSet precutcfg( tmp.get<fcllite::PSet>("LEEPreCut") );
  larlite::LEEPreCut  precutalgo;
  precutalgo.configure( precutcfg );
  precutalgo.initialize(); // actually does nothing

  // Setup app output
  bool make_precut_tree = pset.get<bool>("MakePreCutOutputTree",false);
  std::string precut_tfile_outname = "";
  if ( make_precut_tree )
    precut_tfile_outname = pset.get<std::string>("PreCutOutputTFileName","");
  TFile fout( precut_tfile_outname.c_str(), "recreate");
  TTree tout( "dlprecut", "DL PMT Precut Output Tree");
  precutalgo.bindOutputVariablesToTree( &tout );
  
  // Setup Input Data Coordindator
  larlitecv::DataCoordinator dataco;
  dataco.set_filelist( pset.get<std::string>("LArCVInputFilelist"),   "larcv"   );
  dataco.set_filelist( pset.get<std::string>("LArLiteInputFilelist"), "larlite" );
  dataco.configure( cfg_file, "StorageManager", "IOManager", "ApplyPMTPrecuts" );
  dataco.initialize();

  // Get run entries
  int nentries = dataco.get_nentries("larcv");
  std::cout << "NUMBER OF ENTRIES: LARCV=" << nentries << " LARLITE=" << dataco.get_nentries("larlite") << std::endl;

  int user_nentries   = pset.get<int>("NumEntries",-1);
  int user_startentry = pset.get<int>("StartEntry",-1);
  int startentry = 0;
  if ( user_startentry>=0 ) {
    startentry = user_startentry;
  }
  int endentry = nentries;
  if ( user_nentries>=0 ) {
    if ( user_nentries+startentry<nentries )
      endentry = user_nentries+startentry;
  }
  if ( startentry>=nentries ) {
    std::cout << "Starting beyond end of file. Nothing to do." << std::endl;
    return 0;
  }
  std::cout << "Start Entry: " << startentry << std::endl;
  std::cout << "End Entry: " << endentry-1 << std::endl;
  std::cout << "Buckle up!" << std::endl;

  for (int ientry=startentry; ientry<endentry; ientry++) {
    
    dataco.goto_entry(ientry,"larcv");
    int run,subrun,event;
    dataco.get_id(run,subrun,event);
    // ------------------------------------------------------------------------------------------//
    
    std::cout << "=====================================================" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << " Entry " << ientry << " : " << run << " " << subrun << " " << event << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "=====================================================" << std::endl;
    
    // ------------------------------------------------------------------------------------------//
    // Run PreCut
    bool passes_precut = precutalgo.analyze( &(dataco.get_larlite_io()) );
    
    std::cout << "  PMT precut result: " << passes_precut << std::endl;
    if ( make_precut_tree )
      tout.Fill();

    if ( passes_precut ) {
      dataco.save_entry();      
    }
    
  }

  dataco.finalize();
  fout.Write();
  fout.Close();
  
  return 0;
}
