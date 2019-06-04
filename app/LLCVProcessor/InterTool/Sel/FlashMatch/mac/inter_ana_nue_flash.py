import os,sys,gc

if len(sys.argv) != 9:
    print
    print
    print "SUPERA_FILE = str(sys.argv[1])"
    print "PGRAPH_FILE = str(sys.argv[2])"
    print "TRK_FILE    = str(sys.argv[3])"
    print "FLASH_FILE  = str(sys.argv[4])"
    print "CONFIG      = str(sys.argv[5])"
    print "NUM         = str(sys.argv[6]"
    print "IS_MC       = int(sys.argv[7])"
    print "OUT_DIR     = str(sys.argv[8])"
    print
    print
    sys.exit(1)

SUPERA_FILE = str(sys.argv[1])
PGRAPH_FILE = str(sys.argv[2])
TRK_FILE    = str(sys.argv[3])
FLASH_FILE  = str(sys.argv[4])
CONFIG      = str(sys.argv[5])
NUM         = str(sys.argv[6])
IS_MC       = int(str(sys.argv[7]))
OUT_DIR     = str(sys.argv[8])

import ROOT
from larlitecv import larlitecv
from ROOT import llcv

BASE_PATH = os.path.realpath(__file__)
BASE_PATH = os.path.dirname(BASE_PATH)
print "Base path at: ",BASE_PATH
sys.path.insert(0,BASE_PATH)

# processor
proc = llcv.Processor()

# intermodule
imod = llcv.InterModule()

# configure the driver
driver = imod.Driver()

driver.SetOutputFilename("flash_ana_nue_%s.root" % NUM);

selection = llcv.InterSelNueFlashMatch()
driver.AddSelection(selection);

# process
proc.add_llcv_ana(imod)

proc.configure(CONFIG)
proc.dataco().set_outputfile(os.path.join(OUT_DIR, "aho.root"),"larcv")

proc.add_lcv_input_file(SUPERA_FILE)
proc.add_lcv_input_file(PGRAPH_FILE)
proc.add_ll_input_file(FLASH_FILE)
proc.add_ll_input_file(TRK_FILE)

proc.initialize()

proc.batch_process_lcv_reverse()

proc.finalize()

sys.exit(0)
