#include "LarliteFileManager.h"

namespace larlitecv {

  LarliteFileManager::LarliteFileManager( std::string filelist, bool use_cache) 
    : FileManager( filelist, use_cache )
  {
    // ok
  }

  std::string LarliteFileManager::filetype() {
    return "larlite";
  }

  void LarliteFileManager::user_build_index( const std::vector<std::string>& input, std::map< RSE, int >& rse2entry, std::map< int, RSE >& entry2fse ) {
    
  }
			 


}
