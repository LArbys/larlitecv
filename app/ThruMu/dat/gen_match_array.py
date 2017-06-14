import os,sys
try:
   import cPickle as pickle
except:
   import pickle

""" this script takes the pickled matches and dumps a C++ header file """

use_sce = True # for space charge effects

matches = {}
positions = {}

modes = ["tight","loose"]
positiontypes = [ "top_positions", "bottom_positions", "upstream_positions", "downstream_positions" ]

for mode in modes:
   matchtypes = [ "top_matches", "bottom_matches", "upstream_matches", "downstream_matches" ]
   for matchtype in matchtypes:
      if use_sce:
         fin = open("%s_sce_%s.pickle"%(matchtype,mode),'r')
      else:
         fin = open("%s.pickle"%(matchtype),'r')
      matches[(matchtype,mode)] = pickle.load( fin )
      fin.close()
      print matchtype,mode,": ",len(matches[(matchtype,mode)])

   for postype in positiontypes:
      if use_sce:
         fin = open("%s_sce_%s.pickle"%(postype,mode),'r')
      else:
         fin = open("%s.pickle"%(postype),'r')
      positions[(postype,mode)] = pickle.load( fin )
      fin.close()
      print postype,mode,": ",len(positions[(postype,mode)])

header1 = """
#ifndef __BOUNDARY_MATCH_ARRAY__
#define __BOUNDARY_MATCH_ARRAY__
// This file is generated sing gen_match_array.py

namespace larlitecv {
"""
varstrings = []

for mode in modes:

   s = """
   static const int numTopMatches_%s        = %d;
   static const int numBottomMatches_%s     = %d;
   static const int numUpstreamMatches_%s   = %d;
   static const int numDownstreamMatches_%s = %d;
""" % ( mode,len(matches[("top_matches",mode)]),
        mode,len(matches[("bottom_matches",mode)]),
        mode,len(matches[("upstream_matches",mode)]),
        mode,len(matches[("downstream_matches",mode)]) )
   varstrings.append(s)
for mode in modes:

   for matchtype in matchtypes:
      s  = "    static const int %s_%s[%d][3] =  {\n" % (matchtype,mode,len(matches[(matchtype,mode)]))
      for match in matches[(matchtype,mode)]:
         s += "      {%d,%d,%d},\n"%( match[0],match[1],match[2] )
      s += "    };\n"
      varstrings.append(s)

   for postype in positiontypes:
      s  = "    static const float %s_%s[%d][3] =  {\n" % (postype,mode,len(positions[(postype,mode)]))
      for match in positions[(postype,mode)]:
         s += "      {%.2f,%.2f,%.2f},\n"%( match[0],match[1],match[2] )
      s += "    };\n"
      varstrings.append(s)

header2 = """

    class BoundaryMatchArrays {

       public:
         typedef enum { kTight=0, kLoose } MatchMode_t;
         typedef enum { kTop=0, kBottom, kUpstream, kDownstream } Boundary_t;

         BoundaryMatchArrays( MatchMode_t mode);
         virtual ~BoundaryMatchArrays() {};

         int nmatches(Boundary_t b) { 
             if ( mode==kTight ) {
               switch (b) {
                 case kTop:
                   return numTopMatches_tight;
                   break;
                 case kBottom:
                   return numBottomMatches_tight;
                   break;
                 case kUpstream:
                   return numUpstreamMatches_tight;
                   break;
                 case kDownstream:
                   return numDownstreamMatches_tight;
                   break;
                 default:
                   return 0;
                   break;
                };
              }
              else if ( mode==kLoose ) {
                switch (b) {
                 case kTop:
                   return numTopMatches_loose;
                   break;
                 case kBottom:
                   return numBottomMatches_loose;
                   break;
                 case kUpstream:
                   return numUpstreamMatches_loose;
                   break;
                 case kDownstream:
                   return numDownstreamMatches_loose;
                   break;
                 default:
                   return 0;
                   break;
                };
              }
           };

          void getMatch( Boundary_t b, int imatch, int& uwire, int& vwire, int& ywire ) {
             if ( fMode==kTight ) {
               switch (b) {
                 case kTop:
                   uwire =  top_matches_tight[imatch][0];
                   vwire =  top_matches_tight[imatch][1];
                   ywire =  top_matches_tight[imatch][2];
                   break;
                 case kBottom:
                   uwire =  bottom_matches_tight[imatch][0];
                   vwire =  bottom_matches_tight[imatch][1];
                   ywire =  bottom_matches_tight[imatch][2];
                   break;
                 case kUpstream:
                   uwire =  upstream_matches_tight[imatch][0];
                   vwire =  upstream_matches_tight[imatch][1];
                   ywire =  upstream_matches_tight[imatch][2];
                   break;
                 case kDownstream:
                   uwire =  downstream_matches_tight[imatch][0];
                   vwire =  downstream_matches_tight[imatch][1];
                   ywire =  downstream_matches_tight[imatch][2];
                   break;
                 default:
                   uwire = vwire = ywire = -1;
                   break;
                };
              }
              else if ( fMode==kLoose ) {
               switch (b) {
                 case kTop:
                   uwire =  top_matches_loose[imatch][0];
                   vwire =  top_matches_loose[imatch][1];
                   ywire =  top_matches_loose[imatch][2];
                   break;
                 case kBottom:
                   uwire =  bottom_matches_loose[imatch][0];
                   vwire =  bottom_matches_loose[imatch][1];
                   ywire =  bottom_matches_loose[imatch][2];
                   break;
                 case kUpstream:
                   uwire =  upstream_matches_loose[imatch][0];
                   vwire =  upstream_matches_loose[imatch][1];
                   ywire =  upstream_matches_loose[imatch][2];
                   break;
                 case kDownstream:
                   uwire =  downstream_matches_loose[imatch][0];
                   vwire =  downstream_matches_loose[imatch][1];
                   ywire =  downstream_matches_loose[imatch][2];
                   break;
                 default:
                   uwire = vwire = ywire = -1;
                   break;
                };
              }
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


