#ifndef __LLCVERR_H__
#define __LLCVERR_H__

#include <iostream>
#include <exception>

namespace llcv {

  class llcv_err : public std::exception {
    
  public:
    llcv_err(std::string msg="") : std::exception()
    {
      _msg = "\033[93m";
      _msg += msg;
      _msg += "\033[00m";
    }

    virtual ~llcv_err() throw() {}
    virtual const char* what() const throw()
    { return _msg.c_str(); }

  private:
    std::string _msg;
    
  };
}

#endif
