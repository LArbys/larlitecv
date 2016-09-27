#ifndef __DATA_COORDINATOR__
#define __DATA_COORDINATOR__

#include "DataCoordinator.h"
#include <string>
#include <map>
#include <vector>

// larlite
#include "DataFormat/DataFormatTypes.h"
#include "DataFormat/storage_manager.h"

// larcv
#include "DataFormat/IOManager.h"
#include "Base/PSet.h"

namespace larlitecv {

  class FileManager;

  class DataCoordinator {

  public:

    DataCoordinator();
    virtual ~DataCoordinator();

    // get iomans
    larlite::storage_manager& get_larlite_io() { return larlite_io; };
    larcv::IOManager&         get_larcv_io()   { return larcv_io; };
    void configure( std::string cfgfile, 
		    std::string larlite_cfgname, 
		    std::string larcv_cfgname, std::string coord_cfgname="DataCoordinator" );

    
    // add input files
    void add_inputfile( std::string file, std::string ftype );
    
    // set input filelist (text file with file paths)
    void set_filelist( std::string flist, std::string ftype );

    // set output files
    void set_outputfile( std::string filepath, std::string ftype );

    // nentries
    int get_nentries( std::string ftype );

    // navigation
    void goto_entry( int entry, std::string ftype );
    void goto_event( int run, int subrun, int event );

    // set id
    void set_id( int run, int subrun, int event );

    // saving output
    void save_entry();

    // load after specifying files
    void initialize();
    void finalize();
    void close();

  
  protected:
    
    std::vector< std::string > fManagerList;
    std::map< std::string, FileManager* > fManagers;
    bool fInit;

    std::map< std::string, std::vector<std::string> > user_filepaths;
    std::map< std::string, std::string > user_filelists;
    std::map< std::string, std::string > user_outpath;
    void prepfilelists();

    // storage managers
    larlite::storage_manager larlite_io;
    larcv::IOManager         larcv_io;
    std::map< std::string, std::string > user_ioconfig;

    // configs
    std::string cfgfile;
    larcv::PSet larlite_pset;
    larcv::PSet   larcv_pset;
    void do_larlite_config( larlite::storage_manager& ioman, larcv::PSet& pset );
    larlite::data::DataType_t get_enum_fromstring( std::string name );
  };


}

#endif
