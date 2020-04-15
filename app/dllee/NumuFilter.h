#ifndef __LARLITECV_DLLEE_NUMU_FILTER_H__
#define __LARLITECV_DLLEE_NUMU_FILTER_H__

#include <map>
#include <tuple>
#include <string>
#include <vector>

namespace larlitecv {
namespace dllee {

  class NumuFilter {

  protected:

    NumuFilter();
    virtual ~NumuFilter();

  public:

    static void ReturnFilteredDictionary(std::string inputfile,
                                         std::map<std::tuple<int,int,int>,bool>& rse,
                                         std::map<std::tuple<int,int,int,int>,bool>& rsev,
                                         int verbosity );

  };
  
}
}

#endif
