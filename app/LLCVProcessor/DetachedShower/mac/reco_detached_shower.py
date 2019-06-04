import os,sys

if len(sys.argv) != 7:
    print 
    print "PGRAPH_FILE = str(sys.argv[1])"
    print "HIT_FILE    = str(sys.argv[2])"
    print "MCINFO_FILE = str(sys.argv[3])"
    print "OUTPUT_DIR  = str(sys.argv[4])"
    print "REQPDG      = int(sys.argv[5]) -- always 0"
    print "RECLUSTER   = int(sys.argv[6]) -- not supported"
    print 
    sys.exit(1)


PGRAPH_FILE = str(sys.argv[1])
HIT_FILE    = str(sys.argv[2])
MCINFO_FILE = str(sys.argv[3])
OUTPUT_DIR  = str(sys.argv[4])
REQPDG      = int(sys.argv[5])
RECLUSTER   = int(sys.argv[6])

num = int(os.path.basename(PGRAPH_FILE).split(".")[0].split("_")[-1])

import ROOT
from larlitecv import larlitecv
from ROOT import llcv

BASE_PATH = os.path.realpath(__file__)
BASE_PATH = os.path.dirname(BASE_PATH)
sys.path.insert(0,BASE_PATH)

#
# Make processor
#
proc = llcv.Processor()

#
# handshake
#
dlhs = llcv.DLHandshake()
proc.add_llcv_ana(dlhs)
    
#
# run shower reco
#
from showerRecoSecondShower import SecondShowerReco3D
print "...set REQPDG=%d" % REQPDG
ssshr3d = SecondShowerReco3D(REQPDG)
proc.add_ll_ana(ssshr3d)

#
# configure
#
proc.configure(os.path.join(BASE_PATH,"cfg","shower_reco.cfg"))
proc.add_lcv_input_file(PGRAPH_FILE)
proc.add_ll_input_file(HIT_FILE)
if MCINFO_FILE != "INVALID":
    proc.add_ll_input_file(MCINFO_FILE)

SS = os.path.join(OUTPUT_DIR,"shower_reco_detached_out_%d.root" % num)
print "set output file @=%s" % SS

proc.set_output_ll_name(SS)

proc.initialize()

proc.batch_process_lcv()

proc.finalize()

sys.exit(0)
