import ROOT
from larcv import larcv as lcv

def ChannelMax(plane):
    name = "ChannelMaxPlane%d" % plane
    return lcv.ChannelMax(ROOT.std.string(name))

def Combine():
    return lcv.CombineImages("CombineImages")

def Mask():
    return lcv.SegmentMask("ShowerSegment")

def ShowerImage():
    return lcv.MaskImage("ShowerImage")


