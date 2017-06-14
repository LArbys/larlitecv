#include "FoxTrotTrackerAlgoTypes.h"

namespace larlitecv {

  FoxStep::FoxStep( const std::vector<float>& pos, const std::vector<float>& dir )
    : m_pos(pos), m_dir(dir), m_isgood(true) {
  }

  FoxStep::FoxStep( const std::vector<double>& pos, const std::vector<double>& dir ) {
    m_pos.resize(pos.size(),0);
    m_dir.resize(dir.size(),0);
    for (int i=0; i<(int)pos.size(); i++) {
      m_pos[i] = pos[i];
      m_dir[i] = dir[i];
    }
    m_isgood = true;
  }

}
