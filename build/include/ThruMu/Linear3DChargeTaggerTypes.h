#ifndef __LINEAR3DCHARGETAGGER_H__
#define __LINEAR3DCHARGETAGGER_H__

#include <vector>

namespace larlitecv {

  class PointInfo {
  public:
    std::vector<float> xyz;
    std::vector<int> wire_id;
    int row;
    float tick;
    std::vector<int> cols;
    std::vector<bool> planehascharge;
    std::vector<bool> planehasbadch;
    bool goodpoint;
    int planeswithcharge;
    PointInfo() {
      planeswithcharge = 0;			
      goodpoint = false;
      tick = 0.0;
      row = 0;
    };
    virtual ~PointInfo() {};

    bool operator==(const PointInfo& rhs) const {
      if ( row!=rhs.row ) return false;
      if ( cols.size()!=rhs.cols.size() ) return false;
      for ( int i=0; i<(int)cols.size(); i++ ) {
        if ( cols[i]!=rhs.cols[i] ) return false;
      };
      return true;
    };
  };

  class PointInfoList : public std::vector<PointInfo> {
  public:
    PointInfoList() {
      num_pts_w_allcharge = 0;
      num_pts_w_majcharge = 0;
      num_pts_w_allbadch = 0;
      num_pts_w_allempty = 0;
      num_pts_good = 0;
    };
    virtual ~PointInfoList() {};

#ifndef __CINT__
#ifndef __CLING__
    // we must hide this from ROOT's interpretors
    void emplace( PointInfo&& pt );
#endif
#endif

    float fractionHasChargeWith3Planes() { return float(num_pts_w_allcharge)/float(size()); };
    float fractionHasBadChOn3Planes() { return float(num_pts_w_allbadch)/float(size()); };
    float fractionHasNoChargeOn3Planes() { return float(num_pts_w_allempty)/float(size()); };
    float fractionHasChargeOnMajorityOfPlanes() { return float(num_pts_w_majcharge)/float(size()); };
    float fractionGood() { return float(num_pts_good)/float(size()); };

    int num_pts_w_allcharge;
    int num_pts_w_majcharge;
    int num_pts_w_allbadch;
    int num_pts_w_allempty;
    int num_pts_good;
  };

}

#endif
