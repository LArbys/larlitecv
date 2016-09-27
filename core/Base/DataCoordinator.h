#ifndef __DATA_COORDINATOR__
#define __DATA_COORDINATOR__

#include "DataCoordinator.h"
#include <string>
#include <map>
#include <vector>

namespace larlitecv {

  class FileManager;

  class DataCoordinator {

  public:

    DataCoordinator();
    virtual ~DataCoordinator();
    
    // add input files
    void add_inputfile( std::string file, std::string ftype );
    
    // set input filelist (text file with file paths)
    void set_filelist( std::string flist, std::string ftype );

    // set output files
    void set_outputfile( std::string filepath, std::string ftype );

    // nentries
    int get_nentries( std::string ftype ) { return 0; };

    // navigation
    void set_entry( int entry, std::string ftype ) {};
    void set_eventid( int run, int subrun, int event ) {};

    // load after specifying files
    void initialize();
    void finalize();
  
  protected:
    
    std::vector< std::string > fManagerList;
    std::map< std::string, FileManager* > fManagers;
    bool fInit;

    std::map< std::string, std::vector<std::string> > user_filepaths;
    std::map< std::string, std::string > user_filelists;
    std::map< std::string, std::string > user_outpath;
    void prepfilelists();

  };


}

#endif
