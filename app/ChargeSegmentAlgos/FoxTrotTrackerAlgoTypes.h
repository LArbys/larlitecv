#ifndef __FOX_TROT_TRACKER_ALGO_TYPES_H__
#define __FOX_TROT_TRACKER_ALGO_TYPES_H__

#include <vector>

namespace larlitecv {

  class FoxStep {
  public:
    FoxStep() { m_isgood = true; };
    FoxStep( const std::vector<float>& pos3d, const std::vector<float>& dir3d );
    FoxStep( const std::vector<double>& pos3d, const std::vector<double>& dir3d );    
    virtual ~FoxStep() {};

    const std::vector<float> pos() const { return m_pos; };
    const std::vector<float> dir() const { return m_dir; };
    bool isgood() const { return !m_isgood; };

    int numAlternativeSteps() const { return m_other_valid_next_pos.size(); };
    std::vector<float> popAlternativeStep();

  protected:
    // step position and direction
    std::vector<float> m_pos;
    std::vector<float> m_dir;
    bool m_isgood;

    // stored alternate segments in reverse goodness order.
    // allows us to explore paths in a recursive manner
    std::vector< std::vector<float> > m_other_valid_next_pos;

  };

  typedef std::vector< FoxStep > FoxTrack;

}

#endif
