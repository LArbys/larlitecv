import os,sys

if len (sys.argv) != 7:
    print 
    print 
    print "CONFIG_FILE  = str(sys.argv[1])"
    print "SUPERA_FILE  = str(sys.argv[2])"
    print "SSNET_FILE   = str(sys.argv[3])"
    print "PGRAPH_FILE  = str(sys.argv[4])"
    print "TRK_FILE     = str(sys.argv[5])"
    print "OUTPUT_DIR   = str(sys.argv[6])"
    print 
    sys.exit(1)

CONFIG_FILE  = str(sys.argv[1])
SUPERA_FILE  = str(sys.argv[2])
SSNET_FILE   = str(sys.argv[3])
PGRAPH_FILE  = str(sys.argv[4])
TRK_FILE     = str(sys.argv[5])
OUTPUT_DIR   = str(sys.argv[6])

#num = 1
num = int(os.path.basename(TRK_FILE).split(".")[0].split("_")[-1])

import ROOT
from larlitecv import larlitecv
from ROOT import llcv

BASE_PATH = os.path.realpath(__file__)
BASE_PATH = os.path.dirname(BASE_PATH)
sys.path.insert(0,BASE_PATH)

#
# Configure Processor
#
proc = llcv.Processor()


trk_algo = llcv.TrackPGraphMatch()
proc.add_llcv_ana(trk_algo)

proc.configure(CONFIG_FILE)
proc.dataco().set_outputfile(os.path.join(OUTPUT_DIR, "track_pgraph_match_%d.root" % num),"larcv")
proc.add_lcv_input_file(SUPERA_FILE)
proc.add_lcv_input_file(SSNET_FILE)
proc.add_lcv_input_file(PGRAPH_FILE)
proc.add_ll_input_file(TRK_FILE)
proc.initialize()

#
# run over larcv entries
#
proc.batch_process_lcv()

proc.finalize()

sys.exit(0)
