#ifndef __KALMAN_FILTER_H__
#define __KALMAN_FILTER_H__

// for documentation, see: http://kalman.sourceforge.net/doc/index.html

#include <vector>
#include "DataFormat/Image2D.h"

#include "kalman/kfilter.hpp"


namespace larlitecv {

  class KFStopMu : public ::Kalman::EKFilter<float,0> {
    
    // template arguments:
    // float: underlying data type
    // first index of a vector or matrix
    
  public:
    
    KFStopMu( const std::vector<larcv::Image2D>& imgs );
    virtual ~KFStopMu();

    // required classes to be instantiated

    // prediction
    // x_i = f(x_(i-1), u_(i-1), w_(i-1))
    // x: predicted state
    // u: inputs
    // w: noise

    // measurement
    // z_i = h(z_(i-1), v_(i-1))
    // z: measured state
    // v: measured noise

    void makeProcess(); //< make new prediction vector x
    void makeMeasure(); //< make measurement vector z

    void makeA(); //< Jacobian df/dx
    void makeW(); //< Jacobian df/dw
    void makeH(); //< Jacobian dh/dx
    void makeV(); //< Jacobian dh/dv
    void makeQ(); //< Covariance of process noise
    void makeR(); //< Covariance of measurement noise

    int n;

  protected:
    const std::vector<larcv::Image2D>& fImages;
    float fdriftv;
    float ftrigger_time;
    float fmax_step_cm;
    float fpixheight_ticks;
    float fpixwidth;
  };

  typedef KFStopMu::Vector Vector;
  typedef KFStopMu::Matrix Matrix;

}

#endif
