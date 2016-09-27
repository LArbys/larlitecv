#include "LarliteFileManager.h"
#include "Hashlib2plus/hashlibpp.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>
#include <set>
#include <vector>

namespace larlitecv {

  LarliteFileManager::LarliteFileManager( std::string filelist, bool use_cache) 
    : FileManager( filelist, use_cache )
  {
    // ok
  }

  std::string LarliteFileManager::filetype() {
    return "larlite";
  }

  void LarliteFileManager::user_build_index( const std::vector<std::string>& input,
					     std::vector<std::string>& finallist,
					     std::map< RSE, int >& rse2entry, std::map< int, RSE >& entry2rse ) {
    
    std::set<std::string> producers;
    std::set<std::string> datatypes;
    std::set<std::string> treeflavors;
    std::map< std::string, std::vector<std::string> > flavorfiles;
    std::map< std::string, RSElist > file_rselist;
    std::map< RSElist, std::set<std::string> > rse_flavors;

    // in order to build an event index, we need to get for each file
    //   (1) the run, subrun, event number for each file's entry
    //   (2) the trees inside each file.  this collection of trees defines a 'treeflavor'
    // we then need to sort the files in ascending run,subrun,event order.
    // and finally build the total event to RSE index
    // for each event, we need to make sure that the same set of data products are available
    //    so we will use the largest subset of files given that have a consistent set of data products
    for (int ifile=0; ifile<input.size(); ifile++) {
      
      //  get list of keys in the file. this tells us the types of trees
      std::string fpath = input.at(ifile);
      TFile rfile( fpath.c_str(), "OPEN" );
      int nkeys = rfile.GetListOfKeys()->GetEntries();
      bool found_id_tree = false;
      std::set<std::string> trees;

      for (int ikey=0; ikey<nkeys; ikey++) {
	
	std::string keyname = rfile.GetListOfKeys()->At(ikey)->GetName();
	if ( keyname=="larlite_id_tree" ) found_id_tree = true;
	size_t found1 = keyname.find("_");
	size_t found2 = keyname.find("_",found1+1);
	std::string dtype    = keyname.substr(0,found1);
	std::string producer = keyname.substr(found1+1,found2-found1-1 );
	producers.insert( producer );
	datatypes.insert( dtype );
	trees.insert( keyname );

      }
      
      if ( !found_id_tree )
	continue; // skip this file, we won't know how to index it.

      // make a hash out of the name of tree is the file. will be used to define the flavor of this file
      std::string treehashname = ":";
      for ( std::set<std::string>::iterator it=trees.begin(); it!=trees.end(); it++ ) {
	treehashname += (*it)+":";
      }
      hashwrapper *myWrapper = new md5wrapper();
      std::string treehash = myWrapper->getHashFromString( treehashname.c_str() );
      if ( treeflavors.find(treehash)==treeflavors.end() ) {
	flavorfiles.insert( std::pair< std::string, std::vector<std::string> >( treehash, std::vector<std::string>() ) );
      }
      treeflavors.insert(treehash);
      flavorfiles.find(treehash)->second.push_back( fpath );

      // now we want the RSE for each entry of the tree. we use the id tree to get these
      RSElist fileentry_rse;
      TTree* idtree = (TTree*)rfile.Get( "larlite_id_tree" );
      UInt_t run, subrun, event;
      idtree->SetBranchAddress("_run_id",&run);
      idtree->SetBranchAddress("_subrun_id",&subrun);
      idtree->SetBranchAddress("_event_id",&event);

      long bytes = idtree->GetEntry(0);
      long idtree_entry = 0;
      while ( bytes>0 ) {
	RSE entry( run, subrun, event );
	fileentry_rse.emplace_back( entry );
	bytes = idtree->GetEntry( ++idtree_entry );
      }
      file_rselist.insert( std::pair< std::string, RSElist >( fpath, fileentry_rse ) );

      // we associate an event list to a list of file flavors
      bool found_similar_eventlist = false;
      for ( auto &iter : rse_flavors ) {
	if ( iter.first.isequal( fileentry_rse ) ) {
	  found_similar_eventlist = true;
	  iter.second.insert( treehash );
	}
      }
      if ( !found_similar_eventlist ) {
	std::set<std::string> firstfile;
	firstfile.insert(treehash);
	rse_flavors.insert( std::pair< RSElist, std::set<std::string> >( fileentry_rse, std::move(firstfile) ) );
      }

//       std::cout << "File flavor: " << treehash << " number of events: " << fileentry_rse.size() << ": " 
// 		<< fileentry_rse.run() 
// 		<< " " << fileentry_rse.subrun() 
// 		<< " "  << fileentry_rse.event() << std::endl;
	
    }//end of file list loop

    // ok, we now have maps where
    // flavor to list of files
    // RSE list to a list of files
    
    // we now count how many events each flavor-set has
    std::map< std::set<std::string>, int > numevents_per_flavorset;
    for ( auto& iter : rse_flavors ) {
      // have we already seen this flavor set? (if this doesn't work, can hash the flavor sets first
      if ( numevents_per_flavorset.find( iter.second )==numevents_per_flavorset.end() ) {
	numevents_per_flavorset.insert( std::pair< std::set<std::string>, int >(iter.second, 0 ) );
      }
      numevents_per_flavorset.find( iter.second  )->second += (int)iter.first.size();
    }

    // we choose the flavor set with the most events
    int num_in_maxset = -1;
    std::set<std::string> maxset;
    for ( auto& iter: numevents_per_flavorset ) {
      if ( iter.second>num_in_maxset ) {
	num_in_maxset = iter.second;
	maxset = iter.first;
      }
    }

    // now we finally fill what we've been asked to fill
    finallist.clear();
    rse2entry.clear();
    entry2rse.clear();

    // make filelist
    std::set<RSElist> finalrseset;
    for ( auto &flavorset : maxset ) {
      std::vector<std::string>& files = flavorfiles.find( flavorset )->second;
      for ( auto &file : files ) {
	finallist.push_back( file );
	RSElist& rselist = file_rselist.find( file )->second;
	finalrseset.insert( rselist );
      }
    }
    std::vector< RSElist > finalrse_v;
    for ( auto &rselist : finalrseset ) {
      finalrse_v.push_back( rselist );
    }
    std::sort( finalrse_v.begin(), finalrse_v.end() );

    // make rse dictionaries
    int entrynum = 0;
    for ( auto &rselist : finalrse_v ) {
      for ( auto &rse: rselist ) {
	rse2entry.insert( std::pair< RSE, int >( rse, entrynum ) );
	entry2rse.insert( std::pair< int, RSE >( entrynum, rse ) );
	entrynum++;
      }
    }
    
//     std::cout << "Max flavor set has " << numevents_per_flavorset.find(maxset)->second << " entries. "
// 	      << "Consists of " << maxset.size() << " different tree flavors." << std::endl;
//     std::cout << "Index sizes: " << rse2entry.size() << " vs. entries: "<< entrynum << std::endl;
//     std::cout << "Final file list size: " << ffinallist.size() << std::endl;

    
      
  }//end of user_build_index

}
