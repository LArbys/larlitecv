/**
 * \file llcv_base.h
 *
 * \ingroup SP
 * 
 * \brief Class definition file of llcv_base
 *
 * @author Kazu - Nevis 2015 -- thanks again -vic
 */
#ifndef __LLCV_BASE_H__
#define __LLCV_BASE_H__

#include <vector>
#include "llcv_logger.h"

namespace llcv {
    
  class llcv_base {
    
  public:
    
    llcv_base(const std::string logger_name="llcv_base")
      : _logger(nullptr)
    { _logger = &(::llcv::logger::get(logger_name)); }
    
    llcv_base(const llcv_base &original) : _logger(original._logger) {}
    
    virtual ~llcv_base(){};

    inline const llcv::logger& logger() const
    { return *_logger; }
    
    void set_verbosity(::llcv::msg::Level_t level)
    { _logger->set(level); }

    const std::string& name() const
    { return logger().name(); }
    
  private:
    
    llcv::logger *_logger;   //! don't save
    
  };
}
#endif

/** @} */ // end of doxygen group
