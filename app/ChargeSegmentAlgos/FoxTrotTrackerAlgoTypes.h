#ifndef __FOX_TROT_TRACKER_ALGO_TYPES_H__
#define __FOX_TROT_TRACKER_ALGO_TYPES_H__

#include <vector>

namespace larlitecv {

  class FoxStep {
  public:
    FoxStep() { m_isgood = true; };
    FoxStep( const std::vector<float>& pos3d, const std::vector<float>& dir3d );
    virtual ~FoxStep() {};

    const std::vector<float> pos() const { return m_pos; };
    const std::vector<float> dir() const { return m_dir; };
    bool isgood() const { return !m_isgood; };

  protected:
    std::vector<float> m_pos;
    std::vector<float> m_dir;
    bool m_isgood;

  };

  typedef std::vector< FoxStep > FoxTrack;

}

#endif
