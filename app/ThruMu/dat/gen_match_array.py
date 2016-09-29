import os,sys
try:
   import cPickle as pickle
except:
   import pickle

""" this script takes the pickled matches and dumps a C++ header file """

matches = {}
matchtypes = [ "top_matches", "bottom_matches", "upstream_matches", "downstream_matches" ]
for matchtype in matchtypes:
    fin = open("%s.pickle"%(matchtype),'r')
    matches[matchtype] = pickle.load( fin )
    fin.close()
    print matchtype,": ",len(matches[matchtype])

header1 = """
#ifndef __BOUNDARY_MATCH_ARRAY__
#define __BOUNDARY_MATCH_ARRAY__
// This file is generated sing gen_match_array.py

namespace boundaryalgo {

   static const int numTopMatches        = %d;
   static const int numBottomMatches     = %d;
   static const int numUpstreamMatches   = %d;
   static const int numDownstreamMatches = %d;
""" % ( len(matches["top_matches"]),
        len(matches["bottom_matches"]),
        len(matches["upstream_matches"]),
        len(matches["downstream_matches"]) )

varstrings = []
for matchtype in matchtypes:
   s  = "    static const int %s[%d][3] =  {\n" % (matchtype,len(matches[matchtype]))
   for match in matches[matchtype]:
      s += "      {%d,%d,%d},\n"%( match[0],match[1],match[2] )
   s += "    };\n"
   varstrings.append(s)

header2 = """

    class BoundaryMatchArrays {

       public:
         BoundaryMatchArrays() {};
         virtual ~BoundaryMatchArrays() {};

         typedef enum { kTop=0, kBottom, kUpstream, kDownstream } Boundary_t;

         int nmatches(Boundary_t b) { 
             switch (b) {
               case kTop:
                 return numTopMatches;
                 break;
               case kBottom:
                 return numBottomMatches;
                 break;
               case kUpstream:
                 return numUpstreamMatches;
                 break;
               case kDownstream:
                 return numDownstreamMatches;
                 break;
               default:
                 return 0;
                 break;
              };
           };

          void getMatch( Boundary_t b, int imatch, int& uwire, int& vwire, int& ywire ) {
             switch (b) {
               case kTop:
                 uwire =  top_matches[imatch][0];
                 vwire =  top_matches[imatch][1];
                 ywire =  top_matches[imatch][2];
                 break;
               case kBottom:
                 uwire =  bottom_matches[imatch][0];
                 vwire =  bottom_matches[imatch][1];
                 ywire =  bottom_matches[imatch][2];
                 break;
               case kUpstream:
                 uwire =  upstream_matches[imatch][0];
                 vwire =  upstream_matches[imatch][1];
                 ywire =  upstream_matches[imatch][2];
                 break;
               case kDownstream:
                 uwire =  downstream_matches[imatch][0];
                 vwire =  downstream_matches[imatch][1];
                 ywire =  downstream_matches[imatch][2];
                 break;
               default:
                 uwire = vwire = ywire = -1;
                 break;
              };
         };
           
    };//end of class

}
#endif
"""

f = open("BoundaryMatchArrays.h",'w')
print >>f,header1
for s in varstrings:
   print >>f,s
print >>f,header2
f.close()


