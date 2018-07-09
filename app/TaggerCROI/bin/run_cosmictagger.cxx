#include <iostream>
#include <string>

#include "TaggerCROI/CosmicTagger.h"

// This is still in development. But as you can see its trying to make things cleaner
// Basically moving the code in run_tagger.cxx into its own class and trying to
// make things modular so that it can be easily run in a notebook.
// But official tagger is still run_tagger.cxx

int main( int nargs, const char** argv ) {

  std::string cfg = argv[1];
  
  larlitecv::CosmicTagger tagger( cfg );

  for (int n=0; n<tagger.getNumEntries(); n++) {
    tagger.setEntry(n);
    tagger.clearData();
    bool inputok  = tagger.processInputImages();
    bool endptok  = tagger.findBoundaryEnds();
    bool thrumuok = tagger.findThruMuTracks();
    bool stopmuok = tagger.findStopMu();
    bool croiok   = tagger.findCROI();
    bool writeok  = tagger.writeOutput();
  }
  tagger.finalize();
}
