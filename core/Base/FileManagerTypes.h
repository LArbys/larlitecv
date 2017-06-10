#ifndef __FILEMANAGER_TYPES_H__
#define __FILEMANAGER_TYPES_H__

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>

namespace larlitecv {

  class RSE {
  public:
    RSE() {}
    RSE( int _run, int _subrun, int _event, int _subevent=0 ) {
      run      = _run;
      subrun   = _subrun;
      event    = _event;
      subevent = _subevent;
    }
    virtual ~RSE() {};

    int run;
    int subrun;
    int event;
    int subevent;
    
    // make comparison operators
    bool operator<(const RSE &b) const {
      if ( run<b.run ) return true;
      else if ( run>b.run ) return false;
      // run==b.run
      if ( subrun<b.subrun ) return true;
      else if ( subrun>b.subrun) return false;
      //subrun==b.subrun
      if ( event<b.event) return true;
      else if ( event>b.event ) return false;
      //event==b.event
      if ( subevent<b.subevent) return true;
      else if ( subevent>b.subevent ) return false;
      return false;
    };

    bool operator==(const RSE &b) const {
      if ( run==b.run && subrun==b.subrun && event==b.event && subevent==b.subevent ) return true;
      return false;
    };

    friend std::ostream& operator<<(std::ostream &os, RSE const& );
  };
  
  class RSElist : public std::vector< RSE > {
  public:
    RSElist() 
      {};
    virtual ~RSElist() {};
    bool operator<(const RSElist &b) const {
      if ( run()<b.run() ) return true;
      else if ( run()>b.run() ) return false;
      if ( subrun()<b.subrun() ) return true;
      else if ( subrun()>b.subrun()) return false;
      if ( event()<b.event()) return true;
      else if ( event()>b.event() ) return false;
      return false;      
    }
    bool operator==(const RSElist &b) const {
      if ( run()==b.run() 
	   && subrun()==b.subrun()
	   && event()==b.event() ) return true;
      return false;      
    };
    int run()    const { if (size()>0) return at(0).run;    return -1; };
    int subrun() const { if (size()>0) return at(0).subrun; return -1; };
    int event()  const { if (size()>0) return at(0).event;  return -1; };
    bool isequal( const RSElist& b ) const {
      if ( b.size()!=size() ) return false;
      for (int irse=0;irse<(int)size(); irse++) {
        if ( !(b.at(irse)==at(irse)) ) return false;
      } 
      return true;
    };
  };

}

#endif
