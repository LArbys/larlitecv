import os,sys,gc

if len(sys.argv) != 5:
    print
    print
    print "SSNET_FILE = str(sys.argv[1])"
    print "HIT_FILE   = str(sys.argv[2])"
    print "INTER_FILE = str(sys.argv[3])"
    print "OUT_DIR    = str(sys.argv[4])"
    print
    print
    sys.exit(1)

SSNET_FILE = str(sys.argv[1])
HIT_FILE   = str(sys.argv[2])
INTER_FILE = str(sys.argv[3])
OUT_DIR    = str(sys.argv[4])

import ROOT

from larlitecv import larlitecv
from ROOT import llcv

BASE_PATH = os.path.realpath(__file__)
BASE_PATH = os.path.dirname(BASE_PATH)
sys.path.insert(0,BASE_PATH)

TOP_DIR = os.environ['LARLITECV_BASEDIR']
TOP_DIR = os.path.join(TOP_DIR,"app","LLCVProcessor","InterTool")
sys.path.insert(0,os.path.join(TOP_DIR,"Sel"))

from lib.ssnet_modules import attach_ssnet

proc = llcv.Processor()

# attach ssnet
attach_ssnet(proc)

# intermichel
imod = llcv.InterMichel()

# configure the driver
driver = imod.Driver()
FOUT = os.path.join(OUT_DIR,"michel_ana_%d.root")

num = int(SSNET_FILE.split(".")[0].split("_")[-1])
FOUT = FOUT % num

driver.SetOutputFilename(FOUT)

driver.AttachInterFile(INTER_FILE,"vertex_tree")

selection = llcv.SelMichelStudy()
driver.AddSelection(selection);

# process
proc.add_llcv_ana(imod)

proc.configure(os.path.join(BASE_PATH,"inter_michel.cfg"))
proc.dataco().set_outputfile(os.path.join(OUT_DIR, "michel_trash.root"),"larcv")

proc.add_lcv_input_file(SSNET_FILE)
proc.add_ll_input_file(HIT_FILE)

proc.initialize()
proc.batch_process_lcv_reverse()
proc.finalize()

sys.exit(0)
