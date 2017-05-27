#ifndef __CosmicRetagger_h__
#define __CosmicRetagger_h__

/* ================================================
 * CosmicRetagger
 * 
 * This class takes in tagger output and depending on content
 * reconstitues tagger payload quantities and allows rerunning
 * of certain portions.
 *
 * We inherit from CosmicTagger to reuse many of its methods
 * ================================================= */

#include <string>
#include "CosmicTagger.h"

namespace larlitecv {

  class CosmicRetagger : public CosmicTagger {
  public:

    CosmicRetagger( std::string retagger_cfg, std::string taggerout_larcv_file, std::string taggerout_larlite_file );
    virtual ~CosmicRetagger() {};

    virtual bool processInputImages() override; // we change how we load data, using the input of the tagger out

  };
  
}

#endif
