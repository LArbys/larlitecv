import os,sys
import ROOT as rt

supera_dir="/a/data/amsterdam/tmw23/condor/thrumu/mcc7_bnb_cosmic_v00_p00/supera"
thrumu_dir="/a/data/amsterdam/tmw23/condor/thrumu/mcc7_bnb_cosmic_v00_p00/thrumu"

supera_file = {}
supera_entries = {}
flist = os.listdir( supera_dir )
for f in flist:
    if "root" not in f:
        continue
    print f
    froot = rt.TFile(supera_dir+"/"+f)
    segtree = froot.Get("image2d_segment_tree")
    try:
        nentries = segtree.GetEntries()
    except:
        nentries = 0
    if nentries>0:
        supera_file[ int(f.split(".")[-2].split("_")[-1])  ] = supera_dir+"/"+f
        supera_entries[ int(f.split(".")[-2].split("_")[-1]) ] = nentries

print supera_file
supera_ids = supera_file.keys()
supera_ids.sort()

inputlist = open("inputlist.txt",'w')

thrumu_file = {}
flist = os.listdir( thrumu_dir )
for f in flist:
    if "root" not in f or "larcv" not in f:
        continue
    print f
    fid = int(f.split(".")[-2].split("_")[-1])
    if fid not in supera_ids:
        continue

    froot = rt.TFile( thrumu_dir+"/"+f )
    t = froot.Get("partroi_tpc_tree")
    try:
        nentries = t.GetEntries()
    except:
        nentries = -1
    if nentries!=supera_entries[fid]:
        continue

    print >> inputlist,thrumu_dir+"/"+f
    print >> inputlist,supera_file[fid]

inputlist.close()

