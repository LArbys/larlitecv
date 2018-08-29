#ifndef __FILEMANAGER__
#define __FILEMANAGER__

#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include "FileManagerTypes.h"

namespace larlitecv {
 
  class FileManager {
    
  private:
    FileManager();
    
  public:
    FileManager( std::string filelist, bool use_cache=true );
    virtual ~FileManager() {};

    void setFilelist( std::string flist ) { fFilelist = flist; };
    virtual std::string filetype()=0; //< return name of filetype (e.g. larlite, larcv)
    void initialize();
    void getRSE( int entry, int& run, int& subrun, int& event ) const;
    void getEntry( int run, int subrun, int event, int& entry ) const;
    const std::vector<std::string>& get_final_filelist() const { return ffinallist; };
    int nentries() const { return frse2entry.size(); };
    void sortRSE( bool doit ) { m_sort_rse = doit; };
    bool isSorted() { return m_sort_rse; };

  protected:
    
    virtual void user_build_index( const std::vector<std::string>& input,
				   std::vector<std::string>& finalfilelist,
				   std::map< RSE, int >& rse2entry,
				   std::map< int, RSE >& entry2fse ) = 0; ///< pure virtual function that exists frse2entry and fentry2fse are built
    //virtual void user_build_index( const std::vector<std::string>& input ) = 0;
    void parse_filelist( std::vector<std::string>& flist);         ///< parses the filelist
    std::string get_filelisthash(); ///< create md5 hash from filelist contents
    bool cacheExists( std::string hash );
    void load_from_cache( std::string hash );
    void cache_index( std::string hash );
    std::string printset( const std::set< std::string >& myset );

    bool fUseCache;
    bool isParsed;
    bool m_sort_rse;
    std::string fFilelist;
    std::string fFilelistHash;
    
    std::vector< std::string > ffinallist;
    std::map< RSE, int > frse2entry;
    std::map< int, RSE > fentry2rse;

  };


}

#endif
