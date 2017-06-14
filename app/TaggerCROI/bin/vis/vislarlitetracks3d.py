import os,sys
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph import PlotDataItem
import numpy as np
from larlite import larlite

class PyLArLiteDrawTracks3D:

    def __init__(self):
        pass

    def configure(self,pset):
        self.producer = pset.get("producer_name")

    def visualize( self, larlite_io, larcv_io, rawdigit_io ):
        tracks = []
        event_tracks = larlite_io.get_data( larlite.data.kTrack, self.producer )

        ntracks = event_tracks.size()
        for itrack in range(0,ntracks):
            track = event_tracks.at(itrack)
            npts = track.NumberTrajectoryPoints()
            print "track ",itrack," npts=",npts
            x = np.zeros(npts)
            y = np.zeros(npts)
            z = np.zeros(npts)
            for ipt in range(0,npts):
                point3d = track.LocationAtPoint(ipt)
                x[ipt] = point3d.x()-130.0
                y[ipt] = point3d.y()
                z[ipt] = point3d.z()-505.0

            data = np.vstack( [x,y,z] ).transpose()
            plt =  gl.GLLinePlotItem(pos=data, color=(255,255,255,200), width=3, antialias=True)
            plt.uservisname = "track%d_x%d"%(itrack,int(x[0]))
            tracks.append(plt)
        return tracks
    
