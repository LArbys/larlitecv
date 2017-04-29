#include "FoxTrotTrackerAlgoTypes.h"

namespace larlitecv {

  FoxStep::FoxStep( const std::vector<float>& pos, const std::vector<float>& dir )
    : m_pos(pos), m_dir(dir), m_isgood(false) {
    

  }

}
