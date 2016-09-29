import os,sys
import pyqtgraph as pg
from pyqtgraph import PlotDataItem
import numpy as np
from larcv import larcv

class VisAnnodeHits:

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
               
    
    def __init__(self):
        pass

    def configure(self,pset):
        self.producer = pset.get("stage1_annode_producer")

    def visualize( self, larlite_io, larcv_io, rawdigit_io ):
        event_imgs = larcv_io.get_data( larcv.kProductImage2D, self.producer )

        lcv_imgs = event_imgs.Image2DArray()

        cluster_vecs = [] # output container

        for iimg in xrange(0,lcv_imgs.size()):
            lcv_img = lcv_imgs.at(iimg)

            meta = lcv_img.meta()
            plane = meta.plane()
        
            img = larcv.as_ndarray(lcv_img)

            ends = np.argwhere( img>150 )
            x = np.zeros( len(ends) )
            y = np.zeros( len(ends) )
            for ic,end in enumerate(ends):
                color = VisAnnodeHits.COLORS[ ic%VisAnnodeHits.NCOLORS ]
                x[ic] = meta.pos_x(end[0])
                y[ic] = meta.pos_y(end[1])
            plot = PlotDataItem( x=x, y=y, 
                                 pen=None,
                                 symbolBrush=pg.mkBrush(color=color),
                                 symbol='o',symbolPen=pg.mkPen(color=color,width=0.0) )
            cluster_vecs.append(plot)

        print "number of cluster plots: ",len(cluster_vecs)
        return cluster_vecs
            
