import os,sys
import numpy as np
import ROOT as rt
from math import fabs
try:
   import cPickle as pickle
except:
   import pickle

def load_match_list( wire_pixels ):
   matches = {}
   matchtypes = [ "top_matches", "bottom_matches", "upstream_matches", "downstream_matches" ]
   for matchtype in matchtypes:
      fin = open("%s.pickle"%(matchtype),'r')
      matches[matchtype] = pickle.load( fin )
      fin.close()
      print matchtype,": ",len(matches[matchtype])
      
   # CONVERT FACTORS
      
   # fill the tensors
   downsize_factor = 3456/wire_pixels
   if 3456%wire_pixels!=0:
      downsize_factor+=0
   print "downsize factor: ",downsize_factor

   pixel_lists = {} # dict of plane to list of pixel triples
   for matchtype in matchtypes:
      
      pixel_lists[matchtype] = []
      matchlist = matches[matchtype]
      
      for match in matchlist:
         u = match[0]/downsize_factor
         v = match[1]/downsize_factor
         y = match[2]/downsize_factor
         match = (u,v,y)
         if match not in pixel_lists[matchtype]:
            pixel_lists[matchtype].append( match )
   return pixel_lists

if __name__=="__main__":
   matches = load_match_list( 864 )
   for matchtype,matchlist in matches.items():
      print matchtype,": ",len(matchlist)


