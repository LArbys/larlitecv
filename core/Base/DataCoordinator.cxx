#include "DataCoordinator.h"
#include "FileManager.h"
#include "LarcvFileManager.h"
#include "LarliteFileManager.h"
#include <iostream>
#include <fstream>
#include <assert.h>
#include "Base/LArCVBaseUtilFunc.h"

namespace larlitecv {

  DataCoordinator::DataCoordinator() {
    fManagerList.clear();
    fManagers.clear();
    fInit = false;
    fManagerList.push_back("larlite");
    fManagerList.push_back("larcv");
  }

  DataCoordinator::~DataCoordinator() {
    for ( auto &iter : fManagers ) {
      delete iter.second;
      iter.second = nullptr;
    }
  }

  void DataCoordinator::add_inputfile( std::string file, std::string ftype ) {
    if ( user_filepaths.find( ftype )==user_filepaths.end() ) {
      user_filepaths.insert( std::pair<std::string,std::vector<std::string> >( ftype, std::vector<std::string>() ) );
    }
    user_filepaths[ftype].push_back( file );
  }

  void DataCoordinator::set_filelist( std::string flist, std::string ftype ) {
    if ( user_filelists.find(ftype)==user_filelists.end() ) {
      user_filelists.insert( std::pair<std::string,std::string>(ftype,flist) );
    }
    else {
      user_filelists[ftype] = flist;
    }
  }

  void DataCoordinator::set_outputfile( std::string filepath, std::string ftype ) {
    user_outpath.insert( std::pair< std::string, std::string >( ftype, filepath ) );
  }

  void DataCoordinator::initialize() {
    if ( fInit ) {
      std::cout << "Already initialized!" << std::endl;
      return;
    }

    prepfilelists();

    if ( user_filelists.find("larlite")==user_filelists.end() ) {
      std::cout << "larlite filelists has not been prepared" << std::endl;
      return;
    }
    if ( user_filelists.find("larcv")==user_filelists.end() ) {
      std::cout << "larcv filelists has not been prepared" << std::endl;
      return;
    }

    // create manager instances
    LarliteFileManager* flarlite = new LarliteFileManager( user_filelists["larlite"] );
    LarcvFileManager*     flarcv = new LarcvFileManager( user_filelists["larcv"] );

    // pass the filelists to the managers
    fManagers.insert( std::pair< std::string, FileManager* >( "larlite", flarlite ) );
    fManagers.insert( std::pair< std::string, FileManager* >( "larcv",   flarcv ) );

    // this builds the indices, allowing us to sync the processing
    for (auto &iter : fManagers ) {
      iter.second->initialize();
    }

    // now we setup the iomanagers

    // configure them
    larcv_io.configure( larcv_pset );  // we use the configure function
    do_larlite_config( larlite_io, larlite_pset ); // we have to add one for larlite
    
    // larlite iomanager
    for ( auto const &larlitefile : fManagers["larlite"]->get_final_filelist() )
      larlite_io.add_in_filename( larlitefile );
    // larcv iomanager
    for ( auto const &larcvfile : fManagers["larcv"]->get_final_filelist() )
      larcv_io.add_in_file( larcvfile );

    larcv_io.initialize();
    larlite_io.open();

  }

  void DataCoordinator::close() {
    larlite_io.close();
    larcv_io.reset();
  }
  
  void DataCoordinator::finalize() {
//     for (auto &iter : fManagers ) {
//       iter.second->finalize();
//     }
  }
  
  void DataCoordinator::prepfilelists() {
    /// prepare the filelists for the different managers

    for ( auto &ftype : fManagerList ) {
      // has the manager been given a list of files?
      if ( user_filelists.find( ftype )==user_filelists.end() ) {
	// no file list. has files been specified?
	if ( user_filepaths.find( ftype )!=user_filepaths.end() ) {
	  // files specified. we make a filelist from these paths
	  std::string tempfilename = "last_filelist_"+ftype+".txt";
	  std::ofstream infile( tempfilename.c_str() );
	  for ( auto &fpath : user_filepaths[ftype] ) {
	    infile << fpath << '\n';
	  }
	  infile.close();
	  // set the filelist
	  user_filelists.insert( std::pair< std::string, std::string >( ftype, tempfilename ) );
	}
	else {
	  // nothing specified. ok.
	  user_filelists.insert( std::pair< std::string, std::string >( ftype, "" ) );
	}
      }
    }
  }

  void DataCoordinator::configure( std::string cfgfile, 
				   std::string larlite_cfgname, 
				   std::string larcv_cfgname, std::string coord_cfgname ) {
    // we parse the cfgfile twice. once for larlite, the other for larcv
    larcv::PSet pset_coord = larcv::CreatePSetFromFile( cfgfile, coord_cfgname );
    // get the 

    larlite_pset = pset_coord.get<larcv::PSet>( larlite_cfgname );
    larcv_pset   = pset_coord.get<larcv::PSet>( larcv_cfgname );
    
  }

  larlite::data::DataType_t DataCoordinator::get_enum_fromstring( std::string name ) {
    for (int i=0; i<larlite::data::kDATA_TYPE_MAX; i++) {
      if ( larlite::data::kDATA_TREE_NAME[i]==name )
	return (larlite::data::DataType_t)i;
    }
    return larlite::data::kUndefined;
  }

  void DataCoordinator::do_larlite_config( larlite::storage_manager& ioman, larcv::PSet& pset ) {
    // mode
    int iomode = pset.get<int>("IOMode");
    ioman.set_io_mode( (larlite::storage_manager::IOMode_t)iomode );

    // specified read/write datatypes
    std::vector<std::string> readonlyvars  = pset.get<std::vector<std::string> >( "ReadOnlyDataTypes" );
    std::vector<std::string> readonlyname  = pset.get<std::vector<std::string> >( "ReadOnlyProducers" );
    std::vector<std::string> writeonlyvars = pset.get<std::vector<std::string> >( "WriteOnlyDataTypes" );
    std::vector<std::string> writeonlyname  = pset.get<std::vector<std::string> >( "WriteOnlyProducers" );

    if ( readonlyvars.size()!=readonlyname.size() ) {
      std::cout << "ERROR: number of read-only data types and names are not the same." << std::endl;
      assert(false);
    }
    if ( writeonlyvars.size()!=writeonlyname.size() ) {
      std::cout << "ERROR: number of write-only data types and names are not the same." << std::endl;
      assert(false);
    }

    for ( int i=0; i<readonlyvars.size(); i++ ) {
      std::string rdvar  = readonlyvars.at(i);
      larlite::data::DataType_t datat = get_enum_fromstring( rdvar );
      if ( datat!=larlite::data::kUndefined ) ioman.set_data_to_read( datat, readonlyname.at(i) );
    }

    for ( int i=0; i<writeonlyvars.size(); i++ ) {
      std::string rdvar  = writeonlyvars.at(i);
      larlite::data::DataType_t datat = get_enum_fromstring( rdvar );
      if ( datat!=larlite::data::kUndefined ) ioman.set_data_to_write( datat, writeonlyname.at(i) );
    }

  }

  void DataCoordinator::goto_entry( int entry, std::string ftype_driver ) {
    int run, subrun, event, other_entry;
    if ( ftype_driver=="larlite" ) {
      larlite_io.go_to( entry );
      fManagers["larlite"]->getRSE( entry, run, subrun, event );
      fManagers["larcv"]->getEntry( run, subrun, event, other_entry );
      larcv_io.read_entry( other_entry );
    }
    else if ( ftype_driver=="larcv" ) {
      larcv_io.read_entry( entry );
      fManagers["larcv"]->getRSE( entry, run, subrun, event );
      fManagers["larlite"]->getEntry( run, subrun, event, other_entry );
      larlite_io.go_to( other_entry );
    }
    else {
      std::cout << "not a filetype: " << ftype_driver << std::endl;
      assert(false);
    }
  }


  void DataCoordinator::goto_event( int run, int subrun, int event ) {
    int entry;
    fManagers["larlite"]->getEntry( run, subrun, event, entry );
    larlite_io.go_to( entry );
    fManagers["larcv"]->getEntry( run, subrun, event, entry );
    larcv_io.read_entry( entry );
  }

  int DataCoordinator::get_nentries( std::string ftype ) {
    if ( fManagers.find(ftype)==fManagers.end() ) return 0;
    return fManagers[ftype]->nentries();
  }

}
