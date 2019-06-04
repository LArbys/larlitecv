import sys,os
from larcv import larcv


iom = larcv.IOManager()
iom.add_in_file(sys.argv[1])
iom.initialize()

nentries=iom.get_n_entries()

df_m = {}
df_m['run']    = []
df_m['subrun'] = []
df_m['event']  = []

df_m['proton_x'] = []
df_m['proton_y'] = []
df_m['proton_z'] = []

df_m['electron_x'] = []
df_m['electron_y'] = []
df_m['electron_z'] = []

df_m['muon_x'] = []
df_m['muon_y'] = []
df_m['muon_z'] = []


for i in xrange(nentries):
    print "@i=",i
    iom.read_entry(int(i))
    
    ev_roi = iom.get_data(larcv.kProductROI,"segment")

    df_m['run'].append(int(ev_roi.run()))
    df_m['subrun'].append(int(ev_roi.subrun()))
    df_m['event'].append(int(ev_roi.event()))

    proton_seen = False
    electron_seen = False
    muon_seen = False

    for roi in ev_roi.ROIArray():

        if roi.PdgCode() == 2212 and proton_seen == False:
            proton_seen = True
            
            df_m['proton_x'].append(roi.LastStep().X())
            df_m['proton_y'].append(roi.LastStep().Y())
            df_m['proton_z'].append(roi.LastStep().Z())

        if roi.PdgCode() == 11 and electron_seen == False:
            electron_seen = True
            
            df_m['electron_x'].append(roi.LastStep().X())
            df_m['electron_y'].append(roi.LastStep().Y())
            df_m['electron_z'].append(roi.LastStep().Z())

        if roi.PdgCode() == 13 and muon_seen == False:
            muon_seen = True
            
            df_m['muon_x'].append(roi.LastStep().X())
            df_m['muon_y'].append(roi.LastStep().Y())
            df_m['muon_z'].append(roi.LastStep().Z())


    if proton_seen == False:
        df_m['proton_x'].append(float(-1))
        df_m['proton_y'].append(float(-1))
        df_m['proton_z'].append(float(-1))

    if electron_seen == False:
        df_m['electron_x'].append(float(-1))
        df_m['electron_y'].append(float(-1))
        df_m['electron_z'].append(float(-1))
    
    if muon_seen == False:
        df_m['muon_x'].append(float(-1))
        df_m['muon_y'].append(float(-1))
        df_m['muon_z'].append(float(-1))


import pandas as pd
df = pd.DataFrame(df_m)
df.to_pickle(os.path.basename(sys.argv[1]).split(".")[0] + ".pkl")
sys.exit(0)
