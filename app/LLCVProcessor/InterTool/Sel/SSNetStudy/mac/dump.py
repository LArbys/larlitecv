import pandas as pd
import root_numpy as rn

res = pd.DataFrame(rn.root2array("ssnet_ana_1.root"))

print "track"
print res.ssnet_track_frac_vv.iloc[0]
print ""
print "shower"
print res.ssnet_shower_frac_vv.iloc[0]
