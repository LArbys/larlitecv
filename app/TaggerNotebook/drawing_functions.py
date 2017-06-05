import os,sys
import cv2
import ROOT
from larcv import larcv
from larlitecv import larlitecv
import numpy as np

# Conver Image2D into numpy array usable by OpenCV draw functions
def make_img( img2d, minval, maxval, colormap=cv2.COLORMAP_HOT ):
    imgcp = larcv.Image2D(img2d)
    img_np = np.transpose( larcv.as_ndarray( imgcp ), (1,0) )
    img_np[ img_np<minval ] = minval
    img_np -= minval
    img_np[ img_np>maxval ] = maxval
    img_np *= 255.0/(maxval)
    img_cv = cv2.cvtColor( img_np, cv2.COLOR_GRAY2BGR )
    img_cv = cv2.applyColorMap(img_cv.astype(np.uint8), colormap )
    return img_cv

# Overlay location of flashes
def draw_flashes( img2d, input_data ):
    meta = input_data.img_v.front().meta()
    nflash = 0
    for iset in range(input_data.opflashes_v.size()):
        opflash_v = input_data.opflashes_v[iset]
        for iflash in range(opflash_v.size()):
            opflash = opflash_v[iflash]
            anode_tick = 3200.0 + opflash.Time()/0.5 + 30.0
            cathode_tick = anode_tick + 258.0/0.111436/0.5 - 90.0
            if nflash%4==0:
                text_offset = 0
            elif nflash%4==1:
                text_offset = 250
            elif nflash%3==2:
                text_offset = 3456-600
            elif nflash%4==3:
                text_offset = 3456-300
            #print (iset,iflash),": anode=",anode_tick," cathode=",cathode_tick
            if anode_tick>meta.min_y() and anode_tick<meta.max_y():
                anode_row = meta.row(anode_tick)
                img2d = cv2.line( img2d, (0,anode_row), (meta.cols(),anode_row), (255,0,255), 1 )
                img2d = cv2.putText( img2d, "A (%d,%d) %d"%(iset,iflash,anode_tick), (text_offset,anode_row+5), cv2.FONT_HERSHEY_PLAIN, 2.0, (255,255,255), 2 )
            if cathode_tick>meta.min_y() and cathode_tick<meta.max_y():
                cathode_row = meta.row(cathode_tick)
                img2d = cv2.line( img2d, (0,cathode_row), (meta.cols(),cathode_row), (0,255,255), 1 )
                img2d = cv2.putText( img2d, "C (%d,%d) %d"%(iset,iflash,cathode_tick), (text_offset,cathode_row+5), cv2.FONT_HERSHEY_PLAIN, 2.0, (255,255,255), 2 )
            nflash += 1
    return

def draw_boundary_points( img_cv_v, inputdata, spacepoint_v, labelpoints=False ):
    endtype_colors = {0:(0,0,255),   # top
                      1:(255,0,0),   # bot
                      2:(0,255,255), # upstream
                      3:(255,0,255), # downstream
                      4:(255,140,0), # anode
                      5:(85,107,47), # cathode
                      6:(189,183,107) } # image ends
    boundary_symbol = {0:"T",
                       1:"B",
                       2:"U",
                       3:"D",
                       4:"A",
                       5:"C",
                       6:"E"}
    for ipt in range(spacepoint_v.size()):
        sp = spacepoint_v.at(ipt)
        imgcoords = larcv.UBWireTool.getProjectedImagePixel( sp.pos(), inputdata.img_v.front().meta(), 3 )
        color = endtype_colors[ sp.type() ]
        for p in range(3):
          row = imgcoords[0]
          col = imgcoords[p+1]
          #cv2.circle( img_cv_v[p], (sp.at(p).col, sp.at(p).row), 6, color, -1 )
          cv2.circle( img_cv_v[p], (col,row), 6, color, -1 )
          if labelpoints:
            #cv2.putText( img_cv_v[p], "%s%d"%(boundary_symbol[sp.type()],ipt), (sp.at(p).col+8, sp.at(p).row+5), cv2.FONT_HERSHEY_PLAIN, 2.0, (255,255,255),2 )
            cv2.putText( img_cv_v[p], "%s%d"%(boundary_symbol[sp.type()],ipt), (col+8, row+5), cv2.FONT_HERSHEY_PLAIN, 2.0, (255,255,255),2 )
    return
