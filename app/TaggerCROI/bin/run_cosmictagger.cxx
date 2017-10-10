#include <iostream>
#include <string>

#include "TaggerCROI/CosmicTagger.h"

int main( int nargs, const char** argv ) {

  std::string cfg = argv[1];
  
  larlitecv::CosmicTagger tagger( cfg );

  tagger.setEntry(0);

  bool inputok  = tagger.processInputImages();
  bool endptok  = tagger.findBoundaryEnds();
  bool thrumuok = tagger.findThruMuTracks();
  bool stopmuok = tagger.findStopMu();
  bool croiok   = tagger.findCROI();
  bool writeok  = tagger.writeOutput();
  tagger.finalize();
}
