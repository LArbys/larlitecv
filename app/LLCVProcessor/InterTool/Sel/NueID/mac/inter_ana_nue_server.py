import os,sys,gc,argparse

parser = argparse.ArgumentParser("Run shower reco")
parser.add_argument('-c', '--config', required=True,type=str,help="Configuration file")
parser.add_argument('-mc','--ismc',   action='store_true',default=False,help="Indicate files have MC information")
parser.add_argument('-d', '--debug',  action='store_true',default=False,help="Run in debug mode")
parser.add_argument('-id','--run-id', type=str,default="000",help="Run ID")
parser.add_argument('-od','--out-dir',type=str,default="./",help="Output directory")
parser.add_argument('-re','--reco2d', type=str,required=True,help="larlite Reco2D file with hits")
parser.add_argument('input_larcv',nargs='+',help="list of larcv files to load")

args = parser.parse_args(sys.argv[1:])

# SUPERA_FILE = str(sys.argv[1])
# SSNET_FILE  = str(sys.argv[2])
# TAGGER_FILE = str(sys.argv[3])
# VTX_FILE    = str(sys.argv[4])
# RECO2D_FILE = str(sys.argv[5])
# CONFIG_FILE = str(sys.argv[6])
# NUM         = str(sys.argv[7])
# IS_MC       = int(str(sys.argv[8]))
# OUT_DIR     = str(sys.argv[9])

import ROOT
from larlitecv import larlitecv
from ROOT import llcv

BASE_PATH = os.path.realpath(__file__)
BASE_PATH = os.path.dirname(BASE_PATH)
print "Base path at: ",BASE_PATH
sys.path.insert(0,BASE_PATH)

TOP_DIR = os.environ['LARLITECV_BASEDIR']
TOP_DIR = os.path.join(TOP_DIR,"app","LLCVProcessor","InterTool")
sys.path.insert(0,os.path.join(TOP_DIR,"Sel"))

# load processor modules
from lib.ssnet_modules import attach_ssnet
from lib.dead_modules import attach_dead
from lib.writeout_modules import attach_writeout

# proc
proc = llcv.Processor()

# attach SSNet
attach_ssnet(proc)

# attach dead wire maker
attach_dead(proc)

# attach the larlite writeout
attach_writeout(proc)

# intermodule
imod = llcv.InterModule()

# configure the driver
driver = imod.Driver()

driver.SetOutputFilename(os.path.join( args.out_dir,"nueid_ana_%s.root" % (args.run_id) ));

selection = llcv.SelNueID()
selection.SetIsMC( args.ismc )
if args.debug:
    selection.SetDebug(True)

driver.AddSelection(selection);

# process
proc.add_llcv_ana(imod)

proc.configure(args.config)

for f in args.input_larcv:
    proc.add_lcv_input_file(f.strip())

proc.add_ll_input_file(args.reco2d)

proc.set_output_ll_name(os.path.join(args.out_dir,"nueid_ll_out_%s.root" % (args.run_id) ))
proc.dataco().set_outputfile(os.path.join(args.out_dir, "nueid_lcv_out_%s.root" % (args.run_id)),"larcv")

proc.initialize()

proc.batch_process_lcv_reverse()

proc.finalize()

sys.exit(0)

