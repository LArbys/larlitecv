import ROOT
from larlite import larlite as ll
from ROOT import cluster, cmtool

def CopyDL2DLRaw():
    algo = ll.CopyDL2DLRaw()
    algo.set_input_producer("dl")
    algo.set_output_producer("dlraw")
    return algo;

def ClusterRaw():
    
    algo = ll.SimpleClustererMultiVertex()
    algo.setHitProducer("dlshr")
    algo.setVtxProducer("dl")
    algo.setOutClusProducer("raw")
    algo.setRadius(0.4)
    algo.setCellSize(2.0)
    algo.setUseVertex(True)
    algo.setVtxRadius(2.0)
    algo.setVerbose(False)
    algo.setROIRadius(100.)
    
    return algo

def CorrelateDLRaw():

    algo = ll.DLRawClusterAss()
    algo.SetThreshold(10)
    algo.SetDebug(False)

    return algo

def Merge():

    merger_instance = ll.VertexClusterMerger()

    polar = cmtool.CBAlgoPolar()
    polar.SetUsePairWise(True)
    polar.SetBufferAngle(0.0)
    polar.SetVerbose(False)
    polar.SetMergeTillConverge(True)
    
    vtxalign = cmtool.CBAlgoVtxAlign3()
    vtxalign.SetUsePairWise(False)
    vtxalign.SetVerbose(False)
    vtxalign.SetMergeTillConverge(True)
    vtxalign.SetMaxAngleDiff(12)
    vtxalign.SetMinParOAngle(15)
    vtxalign.SetMaxMergeDist(3) # radial length multiplier
    vtxalign.SetMinNHits(10)
    
    merger_instance.GetManager().AddMergeAlgo(polar)
    merger_instance.GetManager().AddMergeAlgo(vtxalign)
    merger_instance.GetManager().MergeTillConverge(False)
    
    merger_instance.SetVertexProducer("dl")
    merger_instance.SetClusterProducer("raw_filtered")

    return merger_instance
