#include "FileManager.h"
#include "Hashlib2plus/hashlibpp.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "TFile.h"
#include "TTree.h"

namespace larlitecv {

  FileManager::FileManager( std::string filelist, bool use_cache ) {
    isParsed = false;
    fFilelist = filelist;
    fUseCache = use_cache;
  }

  void FileManager::initialize() {

    fFilelistHash = get_filelisthash();
    std::cout << "Hash: " << fFilelistHash << std::endl;

    if ( fUseCache && cacheExists(fFilelistHash) ) {
      load_from_cache(fFilelistHash);
    }
    else {
      // we need to build this instance up
      std::vector<std::string> files;
      parse_filelist(files);   ///< get a vector of string with the filelist
      if ( files.size()>0 ) {
	user_build_index(files,ffinallist,frse2entry,fentry2rse); ///< goes to concrete class function to build event index
	cache_index( fFilelistHash, frse2entry, fentry2rse );
      }
    }

  }

  void FileManager::parse_filelist( std::vector<std::string>& flist ) {

    if ( fFilelist=="" )
      return; // no files

    std::ifstream infile( fFilelist.c_str() );
    std::string line;
    while (std::getline(infile, line)) {
      if ( line!="" )
	flist.push_back(line);
    }
  }

  std::string FileManager::get_filelisthash() {
    // we take the filelist, and build a hash. this will provide a label for the event index cache for this filelist

    hashwrapper *myWrapper = new md5wrapper();
    std::string hash = myWrapper->getHashFromFile( fFilelist.c_str() );
    delete myWrapper;
    return hash;
  }
  
  void FileManager::cache_index( std::string hash, std::map< RSE, int >& rse2entry, std::map< int, RSE >& entry2rse ) {
    system("mkdir -p .pylardcache");
    std::string cachefile = ".pylardcache/"+hash+".root";
    TFile rcache( cachefile.c_str(), "recreate" );
    TTree tcache("entry2rse","RSE to entry map");
    int run, subrun, event;
    tcache.Branch("run",&run,"run/I");
    tcache.Branch("subrun",&subrun,"subrun/I");
    tcache.Branch("event",&event,"event/I");
    int entry=0;
    for ( auto &rse : fentry2rse ) {
      run = rse.second.run;
      subrun = rse.second.subrun;
      event = rse.second.event;
      tcache.Fill();
      entry++;
    }
    //tcache.Write();
    rcache.Write();
  }

  void FileManager::load_from_cache( std::string hash ) {
    std::string cachefile = ".pylardcache/"+hash+".root";
    TFile rcache( cachefile.c_str(), "open" );
    TTree* tcache = (TTree*)rcache.Get("entry2rse");
    int run, subrun, event;
    tcache->SetBranchAddress("run",&run);
    tcache->SetBranchAddress("subrun",&subrun);
    tcache->SetBranchAddress("event",&event);
    ULong_t entry=0;
    long bytes = tcache->GetEntry(entry);
    while (bytes>0) {
      frse2entry.insert( std::pair< RSE, int>( RSE(run, subrun, event), (int)entry ) );
      fentry2rse.insert( std::pair< int, RSE>( (int)entry, RSE(run, subrun, event) ) );
      bytes = tcache->GetEntry(entry);
      entry++;
    }
    tcache->Write();
    rcache.Close();
  }

  std::string FileManager::printset( const std::set< std::string >& myset ) {
    std::string yo = "";
    for ( auto &s : myset ) {
      yo += " " + s;
    }
    return yo;
  }

  void FileManager::getRSE( int entry, int& run, int& subrun, int& event ) const {
    run =  subrun = event = 0;
    if ( fentry2rse.find( entry )!=fentry2rse.end() ) {
      auto iter = fentry2rse.find( entry );
      run    = iter->second.run;
      subrun = iter->second.subrun;
      event  = iter->second.event;
    }
  }

  void FileManager::getEntry( int run, int subrun, int event, int& entry ) const {
    entry = 0;
    RSE rse(run,subrun,event);
    auto iter = frse2entry.find( rse );
    if ( iter!=frse2entry.end() ) {
      entry = iter->second;
    }
  }
    
  
}
