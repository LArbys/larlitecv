#ifndef __DATA_COORDINATOR__
#define __DATA_COORDINATOR__

#include "DataCoordinator.h"
#include <string>
#include <map>
#include <vector>

// larlite
#include "DataFormat/DataFormatTypes.h"
#include "DataFormat/storage_manager.h"
#include "DataFormat/data_base.h"

// larcv
#include "DataFormat/EventBase.h"
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

    // get/set id
    void set_id( int run, int subrun, int event );
    int run();
    int subrun();
    int event();
    void get_id( int& run, int& subrun, int& event );

    // saving output
    void save_entry();

    // load after specifying files
    void initialize();
    void finalize();
    void close();

    // larlite get data command
    larlite::event_base* get_data(const larlite::data::DataType_t type, const std::string& name);

    // larcv get data command
    larcv::EventBase* get_data(const larcv::ProductType_t type, const std::string& producer);
    
    // wrapped commands because python can't resolve function
    larlite::event_base* get_larlite_data( const larlite::data::DataType_t type, const std::string& name) { return get_data( type, name); };
    larcv::EventBase*    get_larcv_data( const larcv::ProductType_t type, const std::string& producer ) { return get_data( type, producer ); };
  
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
    std::map< std::string, int > fIOmodes;
    std::string fLastDriver;
    bool larcv_unused;
    bool larlite_unused;

    // configs
    std::string cfgfile;
    larcv::PSet larlite_pset;
    larcv::PSet   larcv_pset;
    void do_larlite_config( larlite::storage_manager& ioman, larcv::PSet& pset );
    larlite::data::DataType_t get_enum_fromstring( std::string name );
  };

}

#endif
