import os,sys
import numpy as np
import ROOT as rt
from math import fabs
try:
   import cPickle as pickle
except:
   import pickle

"""
Script by Taritree.

Generates the different matchs.
"""

def binsearch( zpos, zlist, tolerance, searchdim, wiredict ):
    """
    Binary search to find matching value within tolerance. Finds first one.
    zpos:  float tuple (y,z)
    zlist: list of (id,y,z)
    tolerance: float 
    """
    imin = 0
    imax = len(zlist)
    #print zlist[200],zlist[10]
    i = imax/2
    while i!=imin and i!=imax:
       # match
       dist = np.sqrt( np.power(zpos[0]-zlist[i][1],2) + np.power(zpos[1]-zlist[i][2],2) )
       # print i,zpos,zlist[i],dist,imin,imax # for debug
       if dist<tolerance:
          return zlist[i]
       # no match yet
       if zpos[searchdim]<zlist[i][searchdim+1]:
          imax = int(i)
          i = (i+imin)/2
       else:
          imin = int(i)
          i = (imax+i)/2
    return None

class WireData:
   """ class stores wire start/end data. exposes data in such a way that search by y and z position is easy """
   def __init__(self, geoinfofile, sortby="pos" ):
      fgeo = rt.TFile( geoinfofile )
      self.wireinfo = fgeo.Get("imagedivider/wireInfo")

      # start/end stored in 2D dict: [plane][wireid]
      self.y_start = {0:[],1:[],2:[]}  # y-sorted
      self.z_start = {0:[],1:[],2:[]}  # z-sorted
      self.y_end   = {0:[],1:[],2:[]}  # y-sorted
      self.z_end   = {0:[],1:[],2:[]}  # z-sorted

      self.wiredict = {0:{},1:{},2:{}}

      for entry in range(0,self.wireinfo.GetEntries()):
         self.wireinfo.GetEntry(entry)
         self.y_start[self.wireinfo.plane].append( (self.wireinfo.wireID, self.wireinfo.start[1], self.wireinfo.start[2]) )
         self.z_start[self.wireinfo.plane].append( (self.wireinfo.wireID, self.wireinfo.start[1], self.wireinfo.start[2]) )
         self.y_end[self.wireinfo.plane].append( (self.wireinfo.wireID, self.wireinfo.end[1], self.wireinfo.end[2]) )
         self.z_end[self.wireinfo.plane].append( (self.wireinfo.wireID, self.wireinfo.end[1], self.wireinfo.end[2]) )
         self.wiredict[self.wireinfo.plane][self.wireinfo.wireID] = { "start":(self.wireinfo.start[1],self.wireinfo.start[2]), 
                                                                      "end":(self.wireinfo.end[1],self.wireinfo.end[2]) }
      fgeo.Close()

      _nwires = []

      for p in range(0,3):
         _nwires.append( len(self.y_start[p]) )
      nwires = np.array( _nwires )

      for p in range(0,3):
         self.y_start[p] = sorted( self.y_start[p], key=lambda x: x[1] )
         self.z_start[p] = sorted( self.z_start[p], key=lambda x: x[2] )
         self.y_end[p]   = sorted( self.y_end[p], key=lambda x: x[1] )
         self.z_end[p]   = sorted( self.z_end[p], key=lambda x: x[2] )


def gen_matches( wireinfo, absolute_constraint, absolute_tolerance, dist_tolerance ):
   """
   generates (U,V,Y) triples that indicate a border crossing

   strategy: loop through u,v,y by a sorted position list and find matches

   input
   -----
   wireinfo: instance of WireData
   absolute_constraint: (y,z) position that the matches must be close to. Provide None if no match required
   absolute_tolerance:  distance from absolute_constraint point
   dist_tolerance: distance the ends can be from one another

   output
   ------
   list of (u,v,y) tuples
   """
   
   matchlist = []

   if absolute_constraint[0] is not None and absolute_constraint[1] is None:
      # a y-value constraint (top and bottom matches, as y is fixed)
      constraint_dim = "Y"
      condim = 1
      plane = 2 # plane that drives search
      # z-sorted list provides ends to match
      swireconstraint = wireinfo.z_start
      ewireconstraint = wireinfo.z_end
      constraint = absolute_constraint[0]
   elif absolute_constraint[1] is not None and absolute_constraint[0] is None:
      # a z-value constraint (top and bottom matches, as y is fixed)
      constraint_dim = "Z"
      condim = 2
      plane = 0 # plane that drives search
      # y-sorted list provides ends to match
      swireconstraint = wireinfo.y_start
      ewireconstraint = wireinfo.y_end
      constraint = absolute_constraint[1]
   else:
      print "Invalid constraint: ",absolute_constraint
      return
      

   # loop over the end points of the Y wires
   combo = swireconstraint[plane]+ewireconstraint[plane]
   nsearch = 0

   for wiredata in combo: # loop over list (wireid,y,z)
      if np.fabs(wiredata[condim]-constraint)>absolute_tolerance:
         continue
      wid = wiredata[0]
      wire = wireinfo.wiredict[plane][wid]
      if np.fabs(wire["start"][condim-1]-wiredata[condim])<absolute_tolerance:
         pos = wire["start"]
         endtype = "start"
      else:
         pos = wire["end"]
         endtype = "end"
      #print "Constraint ",constraint_dim," planeid=",plane," wireid=",wid," conpos=",pos[condim-1]," pos=",pos," endtype=%s"%(endtype)
      
      # loop through other planes, looking for match
      match = [-1,-1,-1]
      match[plane] = wid
      for p in range(0,3):
        if p==plane:
           continue # do not look for matches from plane that drives the loop
        if constraint_dim=="Z" and p==2:
           continue # for upstream/downstream, don't expect any Y-plane matches. so skip it to save time.

        # pass sorted-positions list to binsearch along with full data dictionary
        if constraint_dim=="Y":
           smatchid = binsearch( pos, swireconstraint[p], dist_tolerance, 1, wireinfo.wiredict[p] )
           ematchid = binsearch( pos, ewireconstraint[p], dist_tolerance, 1, wireinfo.wiredict[p] )
        else:
           smatchid = binsearch( pos, swireconstraint[p], dist_tolerance, 0, wireinfo.wiredict[p] )
           ematchid = binsearch( pos, ewireconstraint[p], dist_tolerance, 0, wireinfo.wiredict[p] )
        #print "  match returned: plane=",p," sid=",smatchid," eid=",ematchid
        if smatchid is not None:
           match[p] = smatchid[0]
        elif ematchid is not None:
           match[p] = ematchid[0]

      if constraint_dim=="Z":
         # in this constraint dimension, we expect no matches in the Z-plane. So we add it in to look at the ends
         if pos[1]<500:
            match[2] = 0
         else:
            match[2] = 3455

      if -1 not in match:
         #print "  match found: ",match,"."
         if match not in matchlist:
            matchlist.append( match )

      nsearch += 1
      #if nsearch>=10:
      #   break

   return matchlist


if __name__ == "__main__":

   wireinfo = WireData( "../geoinfo.root", "pos" )

   match_types = [ ("top_matches",(117,None)), 
                   ("bottom_matches",(-115,None)), 
                   ("upstream_matches",(None,0.25)),
                   ("downstream_matches",(None,1036.75)) ]

   # full loop
   for matchtype,constraint in match_types:
      matches = gen_matches( wireinfo, constraint, 2, 0.30 )
      print "Number of ",matchtype,": ",len(matches)
      fout = open( "%s.pickle"%(matchtype), 'w' )
      pickle.dump( matches, fout )
      fout.close()
              
   # for debug
   #top_matches = gen_matches( wireinfo, (117,None), 2, 0.30 )
   #bot_matches = gen_matches( wireinfo, (-115,None), 2, 0.30 )
   #upstream_matches = gen_matches( wireinfo, (None,0.25), 2, 0.30 )
   #downstream_matches = gen_matches( wireinfo, (None,1036.75), 2, 0.30 )
        


