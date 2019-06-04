import pandas as pd
import root_numpy as rn

res = pd.DataFrame(rn.root2array("michel_ana_1.root"))

print res.iloc[0]
