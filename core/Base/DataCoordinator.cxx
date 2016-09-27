#include "DataCoordinator.h"
#include "FileManager.h"
#include "LarcvFileManager.h"
#include "LarliteFileManager.h"
#include <iostream>
#include <fstream>

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

    // pass the filelists
    fManagers.insert( std::pair< std::string, FileManager* >( "larlite", flarlite ) );
    fManagers.insert( std::pair< std::string, FileManager* >( "larcv",   flarcv ) );

    for (auto &iter : fManagers ) {
      iter.second->initialize();
    }
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
  

}
