#include "StopMuSkeleton.h"

#include <stdexcept>
#include <sstream>

#ifndef __CINT__
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif
#endif

#include "DataFormat/Image2D.h"


namespace larlitecv {
    
  larcv::Image2D StopMuSkeleton::skeletonize( const larcv::Image2D& img, const float thresh, const int kernel_size ) {
#ifdef USE_OPENCV
    // convert to cv::Mat
    cv::Mat imgmat = larcv::as_mat_1FC( img );

    // binarize image
    cv::Mat binimg = cv::Mat::zeros( imgmat.rows, imgmat.cols, CV_32FC1 );
    cv::threshold( imgmat, binimg, thresh, 1, cv::THRESH_BINARY );
    std::cout << "threshed" << std::endl;

    // dilate
    cv::Size dilate_kernel_size(kernel_size,kernel_size);
    cv::Mat dilate_kernel = cv::getStructuringElement( cv::MORPH_ELLIPSE, dilate_kernel_size );
    cv::dilate( binimg, binimg, dilate_kernel );

    // skeletonization element
    cv::Size skel_kernel_size(kernel_size,kernel_size);
    cv::Mat skel_kernel = cv::getStructuringElement( cv::MORPH_CROSS, skel_kernel_size );

    // skeletonization img
    cv::Mat skelimg = cv::Mat::zeros( imgmat.rows, imgmat.cols, imgmat.type() );
    int imgsize = imgmat.rows*imgmat.cols;

    bool done = false;
    int iteration = 0;
    while ( !done ) {

      cv::Mat eroded = cv::Mat::zeros( imgmat.rows, imgmat.cols, imgmat.type() );
      cv::Mat temp = cv::Mat::zeros( imgmat.rows, imgmat.cols, imgmat.type() );
      cv::erode( binimg, eroded, skel_kernel );
      cv::dilate( eroded, temp, skel_kernel );
      cv::subtract( binimg, temp, temp );
      cv::bitwise_or( skelimg, temp, skelimg );
      binimg = eroded;
      int img_zeros = imgsize - cv::countNonZero(binimg);

      if ( img_zeros==imgsize ) done = true;
      iteration += 1;
      if ( iteration>100 ) break;
    }
      
    std::cout << "Skeleton made after " << iteration << " iterations." << std::endl;
      
    return larcv::mat_to_image2d( skelimg, img.meta() );
#else

    std::stringstream ss;
    ss << __FILE__ << ":" << __LINE__ << "  use of " << __PRETTY_FUNCTION__ << " requires OpenCV" << std::endl;
    throw std::runtime_error(ss.str());
#endif
    
    larcv::Image2D a(img.meta());
    return a;
  }


}
