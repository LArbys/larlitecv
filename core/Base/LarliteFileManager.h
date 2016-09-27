#ifndef __LARLITEFILEMANAGER__
#define __LARLITEFILEMANAGER__

#include "FileManager.h"
#include <string>
#include <vector>
#include <set>

namespace larlitecv {
  
  class LarliteFileManager : public FileManager {

  public:

    LarliteFileManager( std::string fman, bool use_cache=true );
    virtual ~LarliteFileManager() {};
    
    virtual std::string filetype();
    std::vector<std::string> const& getfilelist() { return ffinallist; };

  protected:
    virtual void user_build_index( const std::vector<std::string>& input,
				   std::vector<std::string>& finalfilelist,
				   std::map< RSE, int >& rse2entry,
				   std::map< int, RSE >& entry2fse );

    std::vector<std::string> ffinallist;
  };
}

#endif
