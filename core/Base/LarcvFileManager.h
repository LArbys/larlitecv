#ifndef __LARCVFILEMANAGER__
#define __LARCVFILEMANAGER__

#include "FileManager.h"
#include <string>
#include <vector>
#include <set>

namespace larlitecv {
  
  class LarcvFileManager : public FileManager {

  public:

    LarcvFileManager( std::string fman, bool use_cache=true );
    virtual ~LarcvFileManager() {};
    
    virtual std::string filetype();
    std::vector<std::string> const& getfilelist() { return ffinallist; };

  protected:
    virtual void user_build_index( const std::vector<std::string>& input,
				   std::map< RSE, int >& rse2entry,
				   std::map< int, RSE >& entry2fse );

    std::vector<std::string> ffinallist;

  };
}

#endif
