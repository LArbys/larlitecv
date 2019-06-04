import sys
import pandas as pd
import root_numpy as rn

res = pd.DataFrame(rn.root2array(str(sys.argv[1])))

idx = int(sys.argv[2])
print res.iloc[idx].to_string()
