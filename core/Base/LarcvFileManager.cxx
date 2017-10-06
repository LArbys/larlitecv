#include "LarcvFileManager.h"
#include "Hashlib2plus/hashlibpp.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <assert.h>
#include "DataFormat/EventBase.h"
#include "DataFormat/EventROI.h"
#include <algorithm>

namespace larlitecv {

  LarcvFileManager::LarcvFileManager( std::string filelist, bool use_cache) 
    : FileManager( filelist, use_cache )
  {
    // ok
  }

  std::string LarcvFileManager::filetype() {
    return "larcv";
  }

  void LarcvFileManager::user_build_index( const std::vector<std::string>& input, 
					   std::vector<std::string>& finallist,
					   std::map< RSE, int >& rse2entry, std::map< int, RSE >& entry2rse ) {
    std::set<std::string> producers;
    std::set<std::string> datatypes;
    std::set<std::string> treeflavors;
    std::map< std::string, std::vector<std::string> > flavorfiles;
    std::map< std::string, RSElist > file_rselist;
    std::map< RSElist, std::set<std::string> > rse_flavors;
    std::map< RSElist, std::vector<std::string> > rse_filelist; 

    // in order to build an event index, we need to get for each file
    //   (1) the run, subrun, event number for each file's entry
    //   (2) the trees inside each file.  this collection of trees defines a 'treeflavor'
    // we then need to sort the files in ascending run,subrun,event order.
    // and finally build the total event to RSE index
    // for each event, we need to make sure that the same set of data products are available
    //    so we will use the largest subset of files given that have a consistent set of data products
    for (size_t ifile=0; ifile<input.size(); ifile++) {
      
      //  get list of keys in the file. this tells us the types of trees
      std::string fpath = input.at(ifile);
      TFile rfile( fpath.c_str(), "OPEN" );
      int nkeys = rfile.GetListOfKeys()->GetEntries();
      bool found_id_tree = false;
      std::string idtreename = "";
      std::string idtreeproducer = "";
      std::string idtreetype = "";
      std::set<std::string> trees;
      
      for (int ikey=0; ikey<nkeys; ikey++) {
	
        std::string keyname = rfile.GetListOfKeys()->At(ikey)->GetName();
      	size_t foundlast = keyname.find_last_of("_");
        std::string tail = keyname.substr(foundlast+1,std::string::npos);
	       if ( tail!="tree") 
	         continue;
         size_t found1 = keyname.find("_");
         std::string dtype    = keyname.substr(0,found1);
         std::string producer = keyname.substr(found1+1,foundlast-found1-1 );
         if ( (!found_id_tree && dtype=="image2d") or 
	      (!found_id_tree && dtype=="partroi") or
	      (!found_id_tree && dtype=="pgraph" ) ) {
          found_id_tree = true;
          idtreename = keyname;
          idtreetype = dtype;
          idtreeproducer = producer;
	  std::cout << "set idtreeproducer: " << idtreeproducer << " dtype=" << idtreetype << std::endl;
        }
        producers.insert( producer );
        datatypes.insert( dtype );
        trees.insert( keyname );
      }
      
      if ( !found_id_tree ) {
        continue; // skip this file, we won't know how to index it.
      }

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
      delete myWrapper;
      treeflavors.insert(treehash);
      flavorfiles.find(treehash)->second.push_back( fpath );
      
      // now we want the RSE for each entry of the tree. we use the id tree to get these
      RSElist fileentry_rse;
      TTree* idtree = (TTree*)rfile.Get( idtreename.c_str() );
      ULong_t run, subrun, event;
      //larcv::EventROI* ev_roi = nullptr;
      //idtree->SetBranchAddress("_run",&run);
      //idtree->SetBranchAddress("_subrun",&subrun);
      //idtree->SetBranchAddress("_event",&event);
      larcv::EventBase* product_ptr = nullptr; 
      if ( idtreetype=="image2d" )
        product_ptr = (larcv::EventBase*)(larcv::DataProductFactory::get().create(larcv::kProductImage2D,idtreeproducer));
      else if ( idtreetype=="partroi" ) 
        product_ptr = (larcv::EventBase*)(larcv::DataProductFactory::get().create(larcv::kProductROI,idtreeproducer));
      else if ( idtreetype=="pgraph" ) 
        product_ptr = (larcv::EventBase*)(larcv::DataProductFactory::get().create(larcv::kProductPGraph,idtreeproducer));
      else {
        throw std::runtime_error( "could not find a LArCV tree to build an index with." );
      }

      std::string brname = idtreetype + "_" + idtreeproducer + "_branch";
      idtree->SetBranchAddress( brname.c_str(), &(product_ptr) );

      long idtree_entry = 0;
      long bytes = idtree->GetEntry(idtree_entry);

      if ( product_ptr==nullptr || !product_ptr->valid() ) {
        std::string msg = "LarcvFileManageer product_ptr not good. Can't build event index using "+idtreetype+" tree with name="+idtreeproducer;
        throw std::runtime_error(msg);
      }

      std::set<RSE> file_entries;
      std::map<RSE,int> duplicate_counter;
      while ( bytes>0 ) {
        run = product_ptr->run();
        subrun = product_ptr->subrun();
        event = product_ptr->event();
        RSE entry( (int)run, (int)subrun, (int)event );
	if ( file_entries.find(entry)!=file_entries.end() ) {
	  // duplicate RSE!
	  RSE duplicate((int)run,(int)subrun,(int)event);
	  if ( duplicate_counter.find(duplicate)==duplicate_counter.end() )
	    duplicate_counter.insert( std::pair<RSE,int>(duplicate,0) );
	  // add subevent number
	  entry.subevent = ++duplicate_counter[duplicate];
	}
	file_entries.insert(entry);
        fileentry_rse.emplace_back( std::move(entry) );
        idtree_entry++;
        bytes = idtree->GetEntry( idtree_entry );
      }
      //std::cout << "last idtree entry: " << idtree_entry << " " << fileentry_rse.size() << std::endl;
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

      // associate rselist to filelist
      auto iter_rse2flist = rse_filelist.find(fileentry_rse);
      if ( iter_rse2flist==rse_filelist.end() ) {
        // make a new filelist
        std::vector<std::string> newflist;
        newflist.push_back( fpath );
        rse_filelist.insert( std::pair< RSElist, std::vector<std::string> >( fileentry_rse, std::move(newflist) ) );
      }
      else {
        iter_rse2flist->second.push_back(fpath);
      }

      std::cout << "File " << fpath << " file-flavor: " << treehash << " number of events: " << fileentry_rse.size() << ": " 
       		<< fileentry_rse.run() 
       		<< " " << fileentry_rse.subrun() 
       		<< " "  << fileentry_rse.event() << std::endl;
      delete product_ptr;
	
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
    std::vector<RSElist> finalrseset;
    for ( auto &flavorset : maxset ) {
      std::vector<std::string>& files = flavorfiles.find( flavorset )->second;
      for ( auto &file : files ) {
        RSElist& rselist = file_rselist.find( file )->second;
        finalrseset.push_back( rselist ); // this use to be a set....
      }
    }
    std::vector< RSElist > finalrse_v;
    for ( auto &rselist : finalrseset ) {
      finalrse_v.push_back( rselist );
    }
    if ( isSorted() )
      sort( finalrse_v.begin(), finalrse_v.end() );

    // make rse dictionaries
    int entrynum = 0;
    for ( auto &rselist : finalrse_v ) {
      for ( auto &rse: rselist ) {
        rse2entry.insert( std::pair< RSE, int >( rse, entrynum ) );
        entry2rse.insert( std::pair< int, RSE >( entrynum, rse ) );
	//std::cout << entrynum << " " << rse << std::endl;	
        entrynum++;
      }
      
      auto iter_rse2flist = rse_filelist.find( rselist );
      for ( auto &fpath : iter_rse2flist->second ) {
        finallist.push_back( fpath ); // we end up resorting
      }
      
    }
    
    // std::cout << "Max flavor set has " << numevents_per_flavorset.find(maxset)->second << " entries. "
    //   << "Consists of " << maxset.size() << " different tree flavors." << std::endl;
    //std::cout << "Index sizes: " << rse2entry.size() << " vs. entries: "<< entrynum << std::endl;
    //std::cout << "Final file list size: " << ffinallist.size() << std::endl;
      
  }//end of user_build_index

}
