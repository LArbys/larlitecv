#include <opencv2/opencv.hpp>

#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"

#include "InterTool_SelNueID/CosmicTag.h"
#include "InterTool_Util/InterImageUtils.h"


#include <iostream>
#include <string>

int main(int argc, char** argv) {
  
  if (argc != 4) {
    std::cout << std::endl;
    std::cout << "argv[1] = input image" << std::endl;
    std::cout << "argv[2] = dead channel image" << std::endl;
    std::cout << "argv[3] = output image name" << std::endl;
    std::cout << std::endl;
    std::exit(1);
  }
  
  std::string argv1(argv[1]);
  std::string argv2(argv[2]);
  std::string argv3(argv[3]);

  auto mat0 = cv::imread(argv1,cv::IMREAD_GRAYSCALE);
  auto mat1 = cv::imread(argv2,cv::IMREAD_GRAYSCALE);
  auto mat2 = mat0.clone();

  llcv::CosmicTag cosmictag;
  
  cosmictag.TagCosmic(mat0,mat1);

  for (const auto& ctor : cosmictag.CosmicContours())
    mat2 = larocv::MaskImage(mat2,ctor,-1,true);

  auto mat3d = llcv::As8UC3(mat2);
  cosmictag.DrawContours(mat3d);

  cv::imwrite(argv3,mat3d);
  
  return 0;
}
