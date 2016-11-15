#include "KFStopMu.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"
#include <cmath>
#include <assert.h>


namespace larlitecv {

  KFStopMu::KFStopMu( const std::vector<larcv::Image2D>& imgs ) : fImages( imgs ) {
    // state vector: (x, vx, y, vy, z, vz )
    // inputs: 0
    // noise variables: sigma in hit position in each plane: (3x2)
    // number of measurements: 4 for time tick + plane positions
    // noise measurements: error in hit position
    setDim(6, 0, 3, 4, 4); // (state vec, inputs, process noise terms, measurements, measurement terms
    //    (X, U, W, Z, V )

    fdriftv = ::larutil::LArProperties::GetME()->DriftVelocity(); // make it cm/usec
    ftrigger_time = 3200;
    fmax_step_cm = 1; // cm
    fpixheight_ticks = imgs.at(0).meta().pixel_height(); // in microseconds
    fpixwidth        = imgs.at(0).meta().pixel_width(); // in microseconds
  }

  KFStopMu::~KFStopMu() {
  }


  void KFStopMu::makeProcess() {
    // Update the state vector: this is the position and velocity of the track
    
    Vector x_(x.size());
    
    // first we calculate the time the particle will travel, in x, equivlanet to the drift time over one time pixel
    float dx = (fpixheight_ticks*0.5)*fdriftv; // [cm] = [pix] * [use/pix] * [cm/usec]
    // then we find out how long this will take the current particle to travel
    float dt = dx/x(1); // [usec] = [cm] / [cm/us]
    float step[3] = {0.};
    float step_size = 0.0;
    for (int i=0; i<3; i++) {
      step[i] = x(2*i+1)*dt;
      step_size += step[i]*step[i];
    }
    step_size = sqrt(step_size);
    // make sure the step size is not too long
    if ( step_size>dx*sqrt(2) ) {
      // too long, shorten step, recalculate 
      std::cout << "shortening step size: " << step_size << " > " << dx*sqrt(2) << std::endl;
      dt *= (step_size)/(dx*sqrt(2));
      step_size = 0.0;
      for (int i=0; i<3; i++) {
	step[i]= x(2*i+1)*dt;
	step_size+= step[i]*step[i];
      }
      step_size =sqrt(step_size);
    }

    // now update the position
    for (int i=0; i<3; i++) {
      x_(i) = x(i) + step[i]; // update position
      x_(2*i+1) = step[i]/dt;    // update "velocity"
    }

    // swap
    x.swap(x_);

  }
  void KFStopMu::makeMeasure() {
    // this are the predicted measurements

    // get the predicted tick position
    float tick = ftrigger_time + x(0)/fdriftv*2.0; // [tick] + [cm]/[cm/usec]*[2 tick/usec]

    // get the predicted wire position
    float pos[3] = { x(0), x(2), x(4) };
    float wid[3] = {0.0};
    for (int p=0; p<3; p++) {
      wid[p] = larutil::Geometry::GetME()->WireCoordinate( pos, p );
    }
    // get pixel positions
    int row = (int)tick/fpixheight_ticks;
    int col[3] = {0};
    for (int p=0; p<3; p++) {
      col[p] = wid[p]/fpixwidth;
    }
    // check if within image
    if ( row<0 || row>=fImages.at(0).meta().rows() ) {
      std::cout << "KFStopMu::makeMeasure(): the predicted row, "  << row << " is outside the image" << std::endl;
      assert(false);
    }
    for (int p=0; p<3; p++) {
      if ( col[p]<0 || col[p]>=fImages.at(0).meta().cols() ) {
	std::cout << "KFStopMu::makeMeasure(): the predicted col in plane=" << p << ", "  << col[p] << " is outside the image" << std::endl;
	assert(false);	
      }
    }
   
    // moving on: set prediction
    z(0) = row;
    z(1) = col[0];
    z(2) = col[1];
    z(3) = col[2];
    
  }

  void KFStopMu::makeA() {
    // calculates df/dx, this is an N x N matrix, where N is the size of the state vector
    int n = x.size();
    for (int i=0; i<n; i++) {
      for (int j=0; j<3; j++) {
	A(i,j) = 0.0;
      }
    }
    A(0,0) = 1.0;
    A(0,1) = 1.0;
    A(1,1) = 1.0;
    A(2,2) = 1.0;
    A(2,3) = 1.0;
    A(3,3) = 1.0;
    A(4,4) = 1.0;
    A(4,5) = 1.0;
    A(5,5) = 1.0;
  }

  void KFStopMu::makeW() {
    // calculates df/dw, this is jacobian of update function and noise process
    int n = x.size();
    for (int i=0; i<n; i++) {
      for (int j=0; j<3; j++) {
	W(i,j) = 0.0;
      }
    }
    // we model noise on top of step direction (coloumb scattering)
    W(1,0) = 1.0; 
    W(3,1) = 1.0;
    W(5,2) = 1.0;
  }

  void KFStopMu::makeH() {
    // calculates dh/dx, this is jacobian of state-to-measurement map and state variables
    int n = x.size();
    int m = 4; // size  of z
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
	H(i,j) = 0.;
      }
    }

    H(0,0) = 2.0/fdriftv; //  dTick/dx
    H(1,1) = -sqrt(3)/0.6; // dU/dy
    H(1,2) = 1.0/0.6;      // dU/dz
    H(2,1) = sqrt(3)/0.6;  // dV/dy
    H(2,2) = 1.0/0.6;      // dV/dz
    H(3,2) = 1.0/0.3;      // dY/dz
    
  }

  void KFStopMu::makeV() {
    // calculate dh/dv, this is jacobian of state-to-measurment map, h, and measurement noise vector
    int n = 4;
    int m = 4;
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
	V(i,j) = 0.0;
      }
    }

    V(0,0) = 1.0;
    V(1,1) = 1.0;
    V(2,2) = 1.0;
    V(3,3) = 1.0;

  }

  void KFStopMu::makeQ() {
    // covariance of state vector noise
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        Q(i,j) = 0.0;
      }
    }
    Q(0,0) = 0.1*0.1;
    Q(1,1) = 0.1*0.1;
    Q(2,2) = 0.1*0.1;
  }

  void KFStopMu::makeR() {
    // covariance of measurement noise
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	R(i,j) = 0.0;
      }
    }
    R(0,0) = 2.0*2.0;
    R(1,1) = 1.0*1.0;
    R(2,2) = 1.0*1.0;
    R(3,3) = 1.0*1.0;
  }

}

