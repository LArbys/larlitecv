#ifndef __FILEMANAGER__
#define __FILEMANAGER__

#include <string>

namespace lartaur {

  class FileManager {
    
  public:
    FileManager( std::string filelist, bool use_cache=true ) {
      isParsed = false;
      fFilelist = filelist;
    };
    virtual ~FileManager() {};

    virtual std::string filetype()=0; //< return name of filetype (e.g. larlite, larcv)

  protected:

    bool isParsed;
    std::string fFilelist;

  };


}

#endif
