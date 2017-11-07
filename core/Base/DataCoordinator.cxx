#include "DataCoordinator.h"
#include "FileManager.h"
#include "LarcvFileManager.h"
#include "LarliteFileManager.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <assert.h>
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/larcv_logger.h"

namespace larlitecv {

  DataCoordinator::DataCoordinator() 
    : larlite_pset("StorageManager","Verbosity: 02 IOMode: 0")
    , larcv_pset("IOManager","Verbosity: 02 IOMode: 0")
  {
    fManagerList.clear();
    fManagers.clear();
    user_filepaths.clear();
    user_filelists.clear();
    user_outpath.clear();
    fInit = false;
    fManagerList.push_back("larlite");
    fManagerList.push_back("larcv");
    larcv_unused = true;
    larlite_unused = true;
  }

  DataCoordinator::~DataCoordinator() {
    for ( auto &iter : fManagers ) {
      delete iter.second;
      iter.second = nullptr;
    }
  }

  void DataCoordinator::add_inputfile( std::string file, std::string ftype ) {
    if ( user_filepaths.find( ftype )==user_filepaths.end() ) {
      user_filepaths.emplace(ftype, std::vector<std::string>());
    }
    auto iter = user_filepaths.find(ftype);
    iter->second.push_back(file);
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
    std::cout << "[DataCoodinator] Initializing" << std::endl;

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
      std::cout << "[DataCoordinator] initializing filemanager for " << iter.first << std::endl;
      iter.second->initialize();
      std::cout << "  " << iter.first << " loading " << iter.second->get_final_filelist().size() << " files." << std::endl;      
    }

    // now we setup the iomanagers

    // configure them
    larcv_io.configure( larcv_pset );  // we use the configure function
    do_larlite_config( larlite_io, larlite_pset ); // we have to add one for larlite
    
    // user may have override output file (for larcv...)
    auto lc_iter = user_outpath.find("larcv");
    if (lc_iter != user_outpath.end()) {
      auto fname = (*lc_iter).second;
      assert(!fname.empty());
      larcv_io.set_out_file(fname);
    }

    // get the iomode for larcv/larlite
    fIOmodes["larcv"]   = (int)larcv_pset.get<int>("IOMode",0);
    fIOmodes["larlite"] = (int)larlite_pset.get<int>("IOMode",0);

    // determine if any of the inputs are unused
    larcv_unused = false;
    larlite_unused = false;
 
    // most obvious tag that is unused: user sets to -1
    if ( fIOmodes["larcv"]   == -1) larcv_unused = true;
    if ( fIOmodes["larlite"] == -1) larlite_unused = true;
     
    // other way is if iomode is read-only, but there are no events provided
    if ( fIOmodes["larcv"]==0  && fManagers["larcv"]->get_final_filelist().empty()  ) larcv_unused = true;
    if ( fIOmodes["larlite"]==0 && fManagers["larlite"]->get_final_filelist().empty()) larlite_unused = true;

    //
    // now load input files
    //
    // larlite iomanager
    if ( !larlite_unused ) {
      for ( auto const &larlitefile : fManagers["larlite"]->get_final_filelist() ) {
	larlite_io.add_in_filename( larlitefile );
      }
      larlite_io.open();
      larlite_io.enable_event_alignment(false);
    }
    // larcv iomanager
    if ( !larcv_unused ) {
      for ( auto const &larcvfile : fManagers["larcv"]->get_final_filelist() ) {
	larcv_io.add_in_file( larcvfile );
      }
      larcv_io.initialize();      
    }

    if ( larlite_unused && larcv_unused ) {
      std::cout << "Both LARCV and LARLITE unused. Must be an error." << std::endl;
      assert(false);
    }
    
    if ( larlite_unused ) std::cout << "[LARLITE unused]" << std::endl;
    if ( larcv_unused )   std::cout << "[LARCV unused]" << std::endl;

  }

  void DataCoordinator::close() {
    larlite_io.close();
    larcv_io.reset();
  }
  
  void DataCoordinator::finalize() {
    if ( !larlite_unused ) larlite_io.close();
    if ( !larcv_unused )   larcv_io.finalize();
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
    std::cout << "Loading pset=" << coord_cfgname << std::endl;
    larcv::PSet pset_head = larcv::CreatePSetFromFile( cfgfile, "cfg" );
    larcv::PSet pset_coord = pset_head.get<larcv::PSet>( coord_cfgname );

    // get the 
    std::cout << "Loading larlite pset=" << larlite_cfgname << std::endl;
    larlite_pset = pset_coord.get<larcv::PSet>( larlite_cfgname );
    std::cout << "Loading larcv pset=" << larcv_cfgname << std::endl;
    larcv_pset   = pset_coord.get<larcv::PSet>( larcv_cfgname );
    
  }

  void DataCoordinator::configure( larcv::PSet& larcv_io_pset, larcv::PSet& larlite_io_pset ) {
    larlite_pset = larlite_io_pset;
    larcv_pset   = larcv_io_pset;
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
    fIOmodes["larlite"] = (int)iomode;

    std::string outfilename = pset.get<std::string>( "OutFileName", "" );
    if ( iomode==1 || iomode==2 ) {
      if (ioman.output_filename().empty()) {
	if ( outfilename.empty()) {
	  std::cout << "Larlite file is set to write mode, but does not have an output file name." << std::endl;
	  assert(false);
	}
	ioman.set_out_filename( outfilename );
	assert(!ioman.output_filename().empty());
      }
    }

    // specified read/write datatypes
    auto readonlyvars  = pset.get<std::vector<std::string> >( "ReadOnlyDataTypes", std::vector<std::string>() );
    auto readonlyname  = pset.get<std::vector<std::string> >( "ReadOnlyProducers", std::vector<std::string>() );
    auto writeonlyvars = pset.get<std::vector<std::string> >( "WriteOnlyDataTypes", std::vector<std::string>() );
    auto writeonlyname = pset.get<std::vector<std::string> >( "WriteOnlyProducers", std::vector<std::string>() );

    if ( readonlyvars.size()!=readonlyname.size() ) {
      std::cout << "ERROR: number of read-only data types and names are not the same." << std::endl;
      assert(false);
    }
    if ( writeonlyvars.size()!=writeonlyname.size() ) {
      std::cout << "ERROR: number of write-only data types and names are not the same." << std::endl;
      assert(false);
    }

    for ( size_t i=0; i<readonlyvars.size(); i++ ) {
      std::string rdvar  = readonlyvars.at(i);
      larlite::data::DataType_t datat = get_enum_fromstring( rdvar );
      if ( datat!=larlite::data::kUndefined ) ioman.set_data_to_read( datat, readonlyname.at(i) );
    }

    for ( size_t i=0; i<writeonlyvars.size(); i++ ) {
      std::string rdvar  = writeonlyvars.at(i);
      larlite::data::DataType_t datat = get_enum_fromstring( rdvar );
      if ( datat!=larlite::data::kUndefined ) ioman.set_data_to_write( datat, writeonlyname.at(i) );
    }

  }

  void DataCoordinator::goto_entry( int entry, std::string ftype_driver ) {
    int run, subrun, event, other_entry;
    fLastDriver = ftype_driver;
    if ( ftype_driver=="larlite" ) {
      if ( larlite_unused ) {
	std::cout << "[larlite unused. goto_entry driven by larlite stopped.]" << std::endl;
	return;
      }
      larlite_io.go_to( entry );
      fManagers["larlite"]->getRSE( entry, run, subrun, event );
      if ( !larcv_unused ) {
	fManagers["larcv"]->getEntry( run, subrun, event, other_entry );
	// std::cout << "given larlite entry=" << entry  << " with "
	// 	  << " rse=(" << run << ", " << subrun << ", " << event << ")"
	// 	  << " corresponds to larcv entry=" << other_entry << std::endl;
	larcv_io.read_entry( other_entry );
      }
    }
    else if ( ftype_driver=="larcv" ) {
      if ( larcv_unused ) {
	std::cout << "[larcv unused. goto_entry driven by larcv stopped.]" << std::endl;
	return;
      }
      larcv_io.read_entry( entry );
      fManagers["larcv"]->getRSE( entry, run, subrun, event );
      if ( !larlite_unused ) {
	fManagers["larlite"]->getEntry( run, subrun, event, other_entry );
	// std::cout << "given larcv entry=" << entry  << " with "
	// 	  << " rse=(" << run << ", " << subrun << ", " << event << ")"
	// 	  << " corresponds to larlite entry=" << other_entry << std::endl;
	larlite_io.go_to( other_entry, false );
      }
    }
    else {
      std::cout << "not a filetype: " << ftype_driver << std::endl;
      assert(false);
    }
    _current_run = run;
    _current_subrun = subrun;
    _current_event = event;
  }


  void DataCoordinator::goto_event( int run, int subrun, int event, std::string ftype_driver ) {
    int entry;
    fLastDriver = ftype_driver;
    if ( !larlite_unused ) {
      fManagers["larlite"]->getEntry( run, subrun, event, entry );
      larlite_io.go_to( entry, false );
      //larlite_io.set_id( run, subrun, event );
    }
    if ( !larcv_unused ) {
      fManagers["larcv"]->getEntry( run, subrun, event, entry );
      larcv_io.read_entry( entry );
      //larcv_io.set_id( run, subrun, event );
    }
    _current_run = run;
    _current_subrun = subrun;
    _current_event = event;
  }

  int DataCoordinator::get_nentries( std::string ftype ) {
    if ( fManagers.find(ftype)==fManagers.end() ) return 0;
    return fManagers[ftype]->nentries();
  }

  void DataCoordinator::save_entry() {

    if ( !larcv_unused )  larcv_io.set_id( _current_run, _current_subrun, _current_event );
    if ( !larlite_unused) larlite_io.set_id( _current_run, _current_subrun, _current_event );

    if ( !larcv_unused )
      larcv_io.save_entry();
    if ( !larlite_unused ) 
      larlite_io.next_event(true);
    // writing done implicitly when event changes for larlite storage_manager
  }

  void DataCoordinator::set_id( int run, int subrun, int event ) {
    if ( !larcv_unused ) larcv_io.set_id( run, subrun, event );
    if ( !larlite_unused )larlite_io.set_id( run, subrun, event );    
    _current_run    = run;
    _current_subrun = subrun;
    _current_event  = event;
  }

  int DataCoordinator::run() {
    if ( fLastDriver=="larlite" ) return larlite_io.run_id();
    else if ( fLastDriver=="larcv" ) return larcv_io.event_id().run();
    return -1;
  }

  int DataCoordinator::subrun() {
    if ( fLastDriver=="larlite" ) return larlite_io.subrun_id();
    else if ( fLastDriver=="larcv" ) return larcv_io.event_id().subrun();
    return -1;
  }

  int DataCoordinator::event() {
    if ( fLastDriver=="larlite" ) return larlite_io.event_id();
    else if ( fLastDriver=="larcv" ) return larcv_io.event_id().event();
    return -1;
  }

  larlite::event_base* DataCoordinator::get_data( const larlite::data::DataType_t type, const std::string& name) {
    return larlite_io.get_data( type, name );
  }

  larcv::EventBase* DataCoordinator::get_data( const larcv::ProductType_t type, const std::string& producer) {
    return larcv_io.get_data( type, producer );
  }

  void DataCoordinator::get_id( int& run, int& subrun, int& event ) {
    if ( fLastDriver==std::string("larlite") ) {
      run    = _current_run;
      subrun = _current_subrun;
      event  = _current_event;
      return;
    }
    else if ( fLastDriver==std::string("larcv") ) {
      run    = _current_run;
      subrun = _current_subrun;
      event  = _current_event;
      return;
    }
    std::stringstream ss;
    ss << __FILE__ << ":" << __LINE__ << " unrecognised driver file type = '" << fLastDriver << "'" << std::endl;
    throw std::runtime_error(ss.str());
  }

}
