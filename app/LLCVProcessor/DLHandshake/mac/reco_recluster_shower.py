import os,sys,argparse

parser = argparse.ArgumentParser("Run shower reco")
parser.add_argument('-c',    '--config',    required=True,type=str,help="Configuration file")
parser.add_argument('-mc',   '--ismc',      action='store_true',default=False,help="Indicate files have MC information")
parser.add_argument('-id',   '--run-id',    type=str,default="000",help="Run ID")
parser.add_argument('-rc',   '--recluster', action='store_true',default=False,help="Recluster flag on.")
parser.add_argument('-od',   '--out-dir',   type=str,default="./",help="Output directory")
parser.add_argument('-llre', '--reco2d',    type=str,required=True,help="larlite RECO2D file with hits")
parser.add_argument('-llmc', '--mcinfo',    type=str,required=False,default='INVALID',help="larlite MCINFO file. Used if MC.")
parser.add_argument('-dqds', '--dqds-file', type=str,required=True,help='dQ/ds table file')
parser.add_argument('input_larcv',nargs='+',help="list of larcv files to load")

args = parser.parse_args(sys.argv[1:])

# if len(sys.argv) != 9:
#     print 
#     print "PGRAPH_FILE = str(sys.argv[1])"
#     print "HIT_FILE    = str(sys.argv[2])"
#     print "MCINFO_FILE = str(sys.argv[3])"
#     print "OUTPUT_DIR  = str(sys.argv[4])"
#     print "REQPDG      = int(sys.argv[5])"
#     print "RECLUSTER   = int(sys.argv[6])"
#     print "IS_MC       = int(sys.argv[7])"
#     print "CONFIG      = str(sys.argv[8])"
#     print 
#     sys.exit(1)


# PGRAPH_FILE = str(sys.argv[1])
# HIT_FILE    = str(sys.argv[2])
# MCINFO_FILE = str(sys.argv[3])
# OUTPUT_DIR  = str(sys.argv[4])
# REQPDG      = int(sys.argv[5])
# RECLUSTER   = int(sys.argv[6])
# IS_MC       = int(sys.argv[7])
# CONFIG      = str(sys.argv[8])

REQPDG = 0

num = int(args.run_id)

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

#
# handshake
#
dlhs = llcv.DLHandshake()
proc.add_llcv_ana(dlhs)


import reclusterElec
if args.recluster:
    #
    # make shower hits
    #
    dlhs = llcv.ShowerHitMaker()
    proc.add_llcv_ana(dlhs)

    #
    # recluster
    #
    # make the raw cluster
    raw_cluster = reclusterElec.ClusterRaw()
    proc.add_ll_ana(raw_cluster)
    
    # correlate raw to dl
    dlraw = reclusterElec.CorrelateDLRaw()
    proc.add_ll_ana(dlraw)

    # recluster seeds
    recluster = reclusterElec.Merge()
    proc.add_ll_ana(recluster)

else:
    #
    # copy dl to dlraw producer
    #
    dlcopy = reclusterElec.CopyDL2DLRaw()
    proc.add_ll_ana(dlcopy)
    
#
# run shower reco
#
from showerRecoDL import DLShowerReco3D
print "...set REQPDG=%d" % REQPDG
if args.ismc:
    dlshr3d = DLShowerReco3D(REQPDG,is_mc=1,shower_table=args.dqds_file)
else:
    dlshr3d = DLShowerReco3D(REQPDG,is_mc=0,shower_table=args.dqds_file)
proc.add_ll_ana(dlshr3d)

#
# configure
#
proc.configure(args.config)

for lcvfile in args.input_larcv:
    proc.add_lcv_input_file( lcvfile )

proc.add_ll_input_file(args.reco2d)

if args.ismc and args.mcinfo != "INVALID":
    proc.add_ll_input_file(args.mcinfo)

proc.set_output_ll_name(os.path.join(args.out_dir,"shower_reco_out_%d.root" % num))
proc.initialize()

#
# run over larcv entries
#
proc.batch_process_lcv()

proc.finalize()

sys.exit(0)
