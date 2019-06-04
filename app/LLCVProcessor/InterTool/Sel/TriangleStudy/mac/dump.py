import pandas as pd
import root_numpy as rn

res = pd.DataFrame(rn.root2array("triangle_ana_cosmic_track_testing_1.root"))

import sys
idx = int(sys.argv[1])
print res.iloc[idx].to_string()
