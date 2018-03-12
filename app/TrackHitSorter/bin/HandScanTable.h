#ifndef __HANDSCANTABLE_H__
#define __HANDSCANTABLE_H__

#include <string>
#include <map>

namespace larlitecv {

  class HandScanTable {
    
  public:
    
    HandScanTable( std::string tabfile );
    ~HandScanTable() {};


    int GetVertexID( int r, int s, int e );
    
  protected:
    
    std::string fTabFile;


  public:
    
    class RSE {
    public:
    RSE(int r, int s, int e) : run(r),subrun(s),event(e) {};
      ~RSE() {};

      int run;
      int subrun;
      int event;

      bool operator<( const RSE& rhs ) const {
	if ( run<rhs.run ) return true;
	else if ( run==rhs.run ) {
	  if ( subrun<rhs.subrun ) return true;
	  else if ( subrun>rhs.subrun ) return false;
	  else if ( subrun==rhs.subrun ) {
	    if ( event<rhs.event ) return true;
	    else return false;
	  }
	}
	return false;
      };
    };


  protected:
    
    std::map<RSE,int> m_map_rse2vertexid;
    
  };
  
}

#endif
