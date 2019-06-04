#include <iostream>

#include "CVUtil/CVUtil.h"
#include "InterTool_Core/InterImageManager.h"
#include <vector>
#include "DataFormat/Image2D.h"

int main() {

  llcv::InterImageManager iim;

  iim.set_verbosity((llcv::msg::Level_t)0);
  
  auto adc_img = larcv::imread_gray("adc_000.jpg");
  std::vector<larcv::Image2D> adc_img2d_v(3);

  iim._vtx_pixel_v.resize(3);

  for(size_t plane=0; plane<3; ++plane) {
    adc_img2d_v[plane] = adc_img;
    iim._vtx_pixel_v[plane] = std::make_pair(367,282);
  }
    
  iim.SetImage(adc_img2d_v,llcv::kADC);

  int cropx = 100;
  int cropy = 100;
  
  auto ocv_img_v = iim.Image<cv::Mat>(llcv::kADC,cropx,cropy);

  cv::imwrite("out1.png",*ocv_img_v.front());

  // llcv::InterImage ii;
  // ii.mat = cv::Mat(10,10,0);
    
  // std::cout << &ii.mat << std::endl;
  // std::cout << ii.get<cv::Mat>() << std::endl;

  return 0;
}
