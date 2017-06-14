#include <iostream>
#include <string>

#include "TaggerCROI/CosmicTagger.h"

int main( int nargs, const char** argv ) {

  std::string cfg = "/home/twongjirad/working/larbys/dllee_unified/larlitecv/app/TaggerNotebook/tagger.cfg";

  larlitecv::CosmicTagger tagger( cfg );

  tagger.setEntry(0);

  bool inputok  = tagger.processInputImages();
  //bool thrumuok = tagger.findThruMu();
  
}
