from __future__ import print_function
import os,sys
import ROOT
from larlite import larlite as fmwk
from recotool import showerreco
from ROOT import protoshower

def getShowerRecoAlgModular(is_mc,shower_dqds_table=None):
    
    if shower_dqds_table is None:
        #'/usr/local/share/dllee_unified/larlite/UserDev/RecoTool/ShowerReco3D/dqds_mc_xyz.txt
        raise ValueError("shower dqds table not specified!")

    alg = showerreco.ShowerRecoAlgModular()
    alg.SetDebug(False)
    alg.SetVerbose(False)
    
    angle3D = showerreco.Angle3DFromVtxQweighted()
    angle3D.setVerbosity(False)

    trunkangle3D = showerreco.TrunkAngle3DFromVtxQweighted()
    trunkangle3D.setVerbosity(False)
    trunkangle3D.setTrunkLength(10)

    startPt = showerreco.VertexIsStartPoint3D()
    startPt.setVerbosity(False)

    energy = showerreco.LinearEnergy()
    energy.setVerbosity(False)

    # implement position-dependent calibration
    energy.CreateResponseMap(20)
    dQdsAVG = 248.
    fin = open(shower_dqds_table,'r')
    #fin = open('/home/vgenty/dqds_mc_xyz.txt','r')
    for line in fin:
        words = line.split()
        x = float(words[0])
        y = float(words[1])
        z = float(words[2])
        q = float(words[3])
        energy.SetResponse(x,y,z,dQdsAVG/q)

    energy.SetElectronLifetime(1e6) # in us DATA value
    energy.SetRecombFactor(0.62)

    if is_mc == 1:
        energy.SetElecGain(200.) # MCC8.3 value
    elif is_mc == 0:
        energy.SetElecGain(240.) # MCC8.3 value
    else:
        print("BAD IS_MC=%s" % str(is_mc))
        sys.exit(1)

    energy.setVerbosity(False)
    energy.SetFillTree(True)

    dqdx = showerreco.dQdxModuleUVY()
    dqdx.setTrunkLength(3.)
    dqdx.setUseTrunkAngle(True)

    length = showerreco.FillLengthUVY()
    length.SetQFraction(0.8)

    alg.AddShowerRecoModule(angle3D)
    alg.AddShowerRecoModule(trunkangle3D)
    alg.AddShowerRecoModule(startPt)
    alg.AddShowerRecoModule(dqdx)
    alg.AddShowerRecoModule(energy)
    alg.AddShowerRecoModule(length)

    alg.PrintModuleList()
    
    return alg

def DefaultShowerReco3D(req_pdg,is_mc,shower_table):
    
    # Create analysis unit
    ana_unit = fmwk.ShowerReco3D()
    
    # require PDG == 11 for PFParticles (?)
    ana_unit.SetRequirePDG11(req_pdg)
    
    # Attach shower reco alg
    sralg = getShowerRecoAlgModular(is_mc,shower_dqds_table=shower_table)
    ana_unit.AddShowerAlgo(sralg)

    return ana_unit

def DLShowerReco3D(req_pdg,is_mc=1,shower_table=None):

    shower_ana_unit=DefaultShowerReco3D(req_pdg,is_mc,shower_table)
    print("Load DefaultShowerReco3D @ ",shower_ana_unit)
    print("... with req_pdg=%d" % int(req_pdg))

    # set ProtoShower Algo to go from data-products to a ProtoShower object
    protoshoweralg = protoshower.ProtoShowerAlgDL()
    protoshoweralg.SetVertexProducer("dlraw")
    shower_ana_unit.GetProtoShowerHelper().setProtoShowerAlg( protoshoweralg )
    
    shower_ana_unit.SetInputProducer("dlraw")
    shower_ana_unit.SetOutputProducer("showerreco")

    return shower_ana_unit


