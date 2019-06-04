import sys, os
import pandas as pd
import root_numpy as rn

FILE = str(sys.argv[1])
IDX  = int(str(sys.argv[2]))

res = pd.DataFrame(rn.root2array(FILE))

print res.iloc[IDX]
