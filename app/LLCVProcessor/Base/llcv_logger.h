/**
 * \file llcv_logger.h
 *
 * \ingroup LArCV
 * 
 * \brief logger utility class definition header file.
 *
 * @author Kazu - Nevis 2015 -- thanks again -vic
 */

#ifndef __LLCV_LOGGER_H__
#define __LLCV_LOGGER_H__

#include <cstdio>
#include <iostream>
#include <map>
#include "LLCVTypes.h"

namespace llcv {

  /**
     \class logger
     \brief Utility class used to show formatted message on the screen.
     A logger class for larcv. Simply shows a formatted colored message on a screen. \n
     A static getter method is provided to create a sharable logger instance (see larcv_base for useage). \n
  */
  class logger{
    
  public:
    
    /// Default constructor
    logger(const std::string& name="no_name")
      : _ostrm(&std::cout)
      , _name(name)
    {}
    
    /// Default destructor
    virtual ~logger(){};
    
  private:
    
    /// ostream
    std::ostream *_ostrm; //!
    
    /// Level
    msg::Level_t _level;
      
    /// Name
    std::string _name;
    
    /// Set of loggers
    static std::map<std::string,llcv::logger> *_logger_m; 

    /// Shared logger
    static llcv::logger* _shared_logger;

    /// Default logger level
    static msg::Level_t _level_default;
    
  public:

    /// Logger's name
    const std::string& name() const { return _name; }

    /// Verbosity level setter
    void set(const msg::Level_t level) { _level = level; }

    /// Verbosity level getter
    msg::Level_t level() const { return _level; }

    /// Comparison operator for static collection of loggers
    inline bool operator<(const logger& rhs) const
    {
      if(_name < rhs.name()) return true;
      if(_name > rhs.name()) return false;
      return false;
    }
    
    /// Getter of a message instance 
    static logger& get(const std::string name)
    {
      if(!_logger_m) _logger_m = new std::map<std::string,llcv::logger>();
      auto iter = _logger_m->find(name);
      if(iter == _logger_m->end()) {
        iter = _logger_m->emplace(name,logger(name)).first;
        iter->second.set(msg::kNORMAL);
      }
      return iter->second;
    };

    static logger& get_shared();

    /// Default logger level getter
    static msg::Level_t default_level() { return _level_default; }
    /// Default logger level setter (only affect future loggers)
    static void default_level(msg::Level_t l) { _level_default = l; }
    /// Force all loggers to change level
    static void force_level(msg::Level_t l)
    {
      default_level(l);
      for(auto& name_logger : *_logger_m) name_logger.second.set(l);
    }
	
    //
    // Verbosity level checker
    //
    inline bool debug   () const { return _level <= msg::kDEBUG;   }
    inline bool info    () const { return _level <= msg::kINFO;    }
    inline bool normal  () const { return _level <= msg::kNORMAL;  }
    inline bool warning () const { return _level <= msg::kWARNING; }
    inline bool error   () const { return _level <= msg::kERROR;   }
    /// Formatted message (simplest)
    std::ostream& send(const msg::Level_t) const;
    /// Formatted message (function name included)
    std::ostream& send(const msg::Level_t level,
		       const std::string& function ) const;
    /// Formatted message (function name + line number)
    std::ostream& send(const msg::Level_t level,
		       const std::string& function,
		       const unsigned int line_num ) const;
    /// Formatted message (function name + line number + file name)
    std::ostream& send(const msg::Level_t level,
		       const std::string& function,
		       const unsigned int line_num,
		       const std::string& file_name) const;
    
  };
}
//
// Compiler macro for saving us from text typing
//
/// Compiler macro for DEBUG message
#define LLCV_DEBUG()    if( logger().debug   () ) logger().send(::llcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
/// Compiler macro for INFO message
#define LLCV_INFO()     if( logger().info    () ) logger().send(::llcv::msg::kINFO,     __FUNCTION__, __LINE__          )
/// Compiler macro for NORMAL message
#define LLCV_NORMAL()   if( logger().normal  () ) logger().send(::llcv::msg::kNORMAL,   __FUNCTION__                    )
/// Compiler macro for WARNING message
#define LLCV_WARNING()  if( logger().warning () ) logger().send(::llcv::msg::kWARNING,  __FUNCTION__                    )
/// Compiler macro for ERROR message
#define LLCV_ERROR()    if( logger().error   () ) logger().send(::llcv::msg::kERROR,    __FUNCTION__, __LINE__          )
/// Compiler macro for CRITICAL message
#define LLCV_CRITICAL()                           logger().send(::llcv::msg::kCRITICAL, __FUNCTION__, __LINE__, __FILE__)

#define LLCV_SDEBUG()    if(llcv::logger::get_shared().debug())    llcv::logger::get_shared().send(::llcv::msg::kDEBUG,    __FUNCTION__,__LINE__,__FILE__)
#define LLCV_SINFO()     if(llcv::logger::get_shared().info())     llcv::logger::get_shared().send(::llcv::msg::kINFO,     __FUNCTION__,__LINE__         )
#define LLCV_SNORMAL()   if(llcv::logger::get_shared().normal())   llcv::logger::get_shared().send(::llcv::msg::kNORMAL,   __FUNCTION__                  )
#define LLCV_SWARNING()  if(llcv::logger::get_shared().warning())  llcv::logger::get_shared().send(::llcv::msg::kWARNING,  __FUNCTION__                  )
#define LLCV_SERROR()    if(llcv::logger::get_shared().error())    llcv::logger::get_shared().send(::llcv::msg::kERROR,    __FUNCTION__,__LINE__         )
#define LLCV_SCRITICAL() llcv::logger::get_shared().send(::llcv::msg::kCRITICAL, __FUNCTION__,__LINE__,__FILE__)
  
/** @} */ // end of doxygen group logger
#endif
