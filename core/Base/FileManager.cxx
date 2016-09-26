#include "FileManager.h"
#include "Hashlib2plus/hashlibpp.h"
#include <fstream>
#include <iostream>

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
      user_build_index(files,frse2entry,fentry2rse); ///< goes to concrete class function to build event index
    }

  }

  void FileManager::parse_filelist( std::vector<std::string>& flist ) {
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
  
}
