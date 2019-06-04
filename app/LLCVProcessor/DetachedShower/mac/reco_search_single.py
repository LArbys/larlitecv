import os,sys,gc

if len(sys.argv) != 4:
    print
    print
    print "SSNET_FILE = str(sys.argv[1])"
    print "VTX_FILE   = str(sys.argv[2])"
    print "OUT_DIR    = str(sys.argv[3])"
    print
    print
    sys.exit(1)


SSNET_FILE = str(sys.argv[1])
VTX_FILE   = str(sys.argv[2])
OUT_DIR    = str(sys.argv[3])

NUM = int(SSNET_FILE.split(".")[0].split("_")[-1])

import ROOT
from larcv import larcv as lcv
from larlite import larlite as ll
from ROOT import llcv

BASE_PATH = os.path.realpath(__file__)
BASE_PATH = os.path.dirname(BASE_PATH)
sys.path.insert(0,BASE_PATH)

from lcv_modules import *

proc = llcv.Processor()

#
# run ssnet
#

# max
for plane in xrange(3):
    res = ChannelMax(plane)
    proc.add_lc_proc(res)

# combine
combine = Combine()
proc.add_lc_proc(combine)

# mask
mask = Mask()
proc.add_lc_proc(mask)

# image
shr_image = ShowerImage()
proc.add_lc_proc(shr_image)

#
# find single detached shower
#

# algo
single_shr = llcv.SearchAlgoSingle()
single_shr.set_verbosity(2)

# driver
locate_shr = llcv.SearchDetached()
single_shr.set_verbosity(2)

locate_shr.SetAlgo(single_shr)

proc.add_llcv_ana(locate_shr)

#
# process
#
proc.configure(os.path.join(BASE_PATH,"cfg","second_shower.cfg"))
proc.dataco().set_outputfile(os.path.join(OUT_DIR, "search_single_%d.root" % NUM),"larcv")
proc.add_lcv_input_file(SSNET_FILE)
proc.add_lcv_input_file(VTX_FILE)
proc.initialize()

proc.batch_process_lcv_reverse()

proc.finalize()

sys.exit(0)

