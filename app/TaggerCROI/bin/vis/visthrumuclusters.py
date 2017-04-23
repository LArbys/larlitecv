import os,sys
import pyqtgraph as pg
from pyqtgraph import PlotDataItem
import numpy as np
from larcv import larcv

class VisThruMuTrackClusters:

    COLORS = [ (76,0,153,125), # purple
               (102,255,102),  # green
               (51,255,255),   # cyan
               (255,128,0),    # orange
               (0,76,153),     # blue
               (204,0,0),      # red
               (0,102,0),      # dark-green
               (102,0,51),     # dark-magenta
               (255,255,0),    # yellow
               ]
    NCOLORS = len(COLORS)
    PLANE_COLORS = {0:(125,0,0,125),
                    1:(0,125,0,125),
                    2:(0,0,125,125)}
    CH_MARKER = { 0:"t",
                  1:"o",
                  2:"s",
                  3:"d" }               
    
    def __init__(self):
        pass

    def configure(self,pset):
        self.img_producer = pset.get("ImageProducer")
        self.track_producer = pset.get("TrackPixelProducer")
        self.plane_selection = pset.get("plane",-1)
        strcolor = pset.get("color").split("[")[-1].split("]")[0]
        self.color = []
        for c in strcolor.split(","):
            self.color.append( int(c) )
            

    def visualize( self, larlite_io, larcv_io, rawdigit_io ):
        track_plots = []
        event_imgs = larcv_io.get_data( larcv.kProductImage2D, self.img_producer )
        if event_imgs.Image2DArray().size()==0:
            print "no images? stop boundary hit vis."
            return
        meta = event_imgs.Image2DArray().at(0).meta()

        event_tracks = larcv_io.get_data( larcv.kProductPixel2D, self.track_producer )
        for iplane in range(0,3):
            if self.plane_selection>=0 and iplane!=self.plane_selection:
                continue

            tracks = event_tracks.Pixel2DClusterArray( iplane )
            ntracks = tracks.size()
            for itrack in range(0,ntracks):
                track = tracks.at(itrack)
                npts = track.size()
                if npts==0:
                    continue
                print " [visthrumuclusters] preparing plot for thrumu track size=",npts
                x = np.zeros( npts )
                y = np.zeros( npts )
                for ipt in range(0,npts):
                    pixel = track.at(ipt)
                    x[ipt] = meta.pos_x( pixel.X() )
                    y[ipt] = meta.pos_y( pixel.Y() )
                plot = PlotDataItem( x=x, y=y, pen=None,symbolPen=tuple(self.color),symbol='s',symbolBrush=None,symbolSize=1,pxMode=False  )
                plot.uservisname = "p%d_cluster_%s_%d"%(iplane,self.track_producer,itrack)
                track_plots.append( plot )

        print "Number of track plots: ",len(track_plots)
        return track_plots
            
