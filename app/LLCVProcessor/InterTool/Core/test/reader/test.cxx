#include <iostream>

#include "InterTool_Core/InterTTreeHandler.h"

int main() {

  std::string fname = "/home/vgenty/sw/larcv/app/LArOpenCVHandle/ana/intermediate_file/test/intermediate_file_342119397.root";
  std::string tname = "vertex_tree";

  llcv::InterTTreeHandler handler(fname,tname);
  handler.Initialize();
    
  // TFile tf(fname.c_str());
  // TTreeReader ttr(tname.c_str(),&tf);
  // TTreeReaderValue<int> a(ttr,"event");
  // for(size_t i=0; i<20; ++i)  {
  //   ttr.SetEntry(i);
  //   ttr.Next();
  //   std::cout << *a << std::endl;
  // }


  
  return 0;
}
