#ifndef __FILEMANAGER__
#define __FILEMANAGER__

#include <string>
#include <vector>
#include <map>

namespace larlitecv {

  class RSE {
  public:
    RSE() {}
    RSE( int _run, int _subrun, int _event ) {
      run = _run;
      subrun = _subrun;
      event = _event;
    }
    virtual ~RSE() {};

    int run;
    int subrun;
    int event;
    
    // make comparison operators
    bool operator<(const RSE &b) const {
      if ( run<b.run ) return true;
      else if ( run>b.run ) return false;
      if ( subrun<b.subrun ) return true;
      else if ( subrun>b.subrun) return false;
      if ( event<b.event) return true;
      else if ( event>b.event ) return false;
      return false;
    };

    bool operator==(const RSE &b) const {
      if ( run==b.run && subrun==b.subrun && event==b.event ) return true;
      return false;
    };
  };
 
  class FileManager {
    
  private:
    FileManager();
    
  public:
    FileManager( std::string filelist, bool use_cache=true );
    virtual ~FileManager() {};

    void setFilelist( std::string flist ) { fFilelist = flist; };
    virtual std::string filetype()=0; //< return name of filetype (e.g. larlite, larcv)
    void initialize();

  protected:
    
    virtual void user_build_index( const std::vector<std::string>& input,
				   std::map< RSE, int >& rse2entry,
				   std::map< int, RSE >& entry2fse ) = 0; ///< pure virtual function that exists frse2entry and fentry2fse are built
    //virtual void user_build_index( const std::vector<std::string>& input ) = 0;
    void parse_filelist( std::vector<std::string>& flist);         ///< parses the filelist
    std::string get_filelisthash(); ///< create md5 hash from filelist contents
    bool cacheExists( std::string hash ) { return false; };
    void load_from_cache( std::string hash ) {};

    bool fUseCache;
    bool isParsed;
    std::string fFilelist;
    std::string fFilelistHash;

    std::map< RSE, int > frse2entry;
    std::map< int, RSE > fentry2rse;

  };


}

#endif
