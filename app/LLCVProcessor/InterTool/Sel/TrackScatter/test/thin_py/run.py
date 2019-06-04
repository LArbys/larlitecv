import pandas as pd
import root_numpy as rn
import numpy as np

nue_df = pd.DataFrame(rn.root2array("../../mac/nue_track_scatter_ana_1.root"))
res = nue_df.query("run==1&subrun==18&event==879")

tid=0
s3d0 = np.vstack((res.shower_x_vv.str[tid].values[0],res.shower_y_vv.str[tid].values[0],res.shower_z_vv.str[tid].values[0])).T

tid=1
s3d1 = np.vstack((res.shower_x_vv.str[tid].values[0],res.shower_y_vv.str[tid].values[0],res.shower_z_vv.str[tid].values[0])).T

data0 = s3d0.astype(np.float32)
data1 = s3d1.astype(np.float32)
data = np.concatenate((data0,data1))

from larlitecv import larlitecv

import ROOT
from ROOT import llcv
from ROOT import larocv

skel = llcv.Skeletonize()

vtx3d_v = ROOT.std.vector("larocv::data::Vertex3D")()
for ix, d in enumerate(data):
    
    vtx3d = larocv.data.Vertex3D()
    vtx3d.x = d[0]
    vtx3d.y = d[1]
    vtx3d.z = d[2]

    vtx3d_v.push_back(vtx3d)


skel.Initialize(vtx3d_v,0.5)
ret = skel.Run()

print ret.size()
