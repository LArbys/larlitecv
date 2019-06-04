import pandas as pd
import root_numpy as rn

df = pd.DataFrame(rn.root2array("track_pgraph_match_1.root",treename="TrackPGraphMatch"))
print df.iloc[0].to_string()
