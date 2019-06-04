import ROOT
from larlite import larlite as fmwk
from recotool import showerreco
from ROOT import protoshower

def getShowerRecoAlgModular():
    
    alg = showerreco.ShowerRecoAlgModular()
    alg.SetDebug(False)
    alg.SetVerbose(False)
    
    angle3D = showerreco.Angle3DFromVtxQweighted()
    angle3D.setVerbosity(False)

    #startPt = showerreco.VertexIsStartPoint3D()
    startPt = showerreco.NearestStartPoint3D()
    startPt.setVerbosity(False)

    energy = showerreco.LinearEnergy()
    energy.setVerbosity(False)

    # implement position-dependent calibration
    energy.CreateResponseMap(20)
    dQdsAVG = 248.
    #fin = open('/usr/local/share/dllee_unified/larlite/UserDev/RecoTool/ShowerReco3D/dqds_mc_xyz.txt','r')
    fin = open('/home/vgenty/dqds_mc_xyz.txt','r')
    for line in fin:
        words = line.split()
        x = float(words[0])
        y = float(words[1])
        z = float(words[2])
        q = float(words[3])
        energy.SetResponse(x,y,z,dQdsAVG/q)

    energy.SetElectronLifetime(1e6) # in us DATA value
    energy.SetRecombFactor(0.62)
    #energy.SetElecGain(243.) # MCC8.0 data
    energy.SetElecGain(200.) # MCC8.0 value
    energy.setVerbosity(False)
    energy.SetFillTree(True)

    dqdx = showerreco.dQdxModuleUVY()
    dqdx.setTrunkLength(3.)

    length = showerreco.FillLengthUVY()

    alg.AddShowerRecoModule(angle3D)
    alg.AddShowerRecoModule(startPt)
    alg.AddShowerRecoModule(dqdx)
    alg.AddShowerRecoModule(energy)
    alg.AddShowerRecoModule(length)

    alg.PrintModuleList()
    
    return alg

def DefaultShowerReco3D(req_pdg):
    
    # Create analysis unit
    ana_unit = fmwk.ShowerReco3D()
    
    # require PDG == 11 for PFParticles (?)
    ana_unit.SetRequirePDG11(req_pdg)
    
    # Attach shower reco alg
    sralg = getShowerRecoAlgModular()
    ana_unit.AddShowerAlgo(sralg)

    return ana_unit

def SecondShowerReco3D(req_pdg):

    shower_ana_unit=DefaultShowerReco3D(req_pdg)
    print "Load DefaultShowerReco3D @ ",shower_ana_unit
    print "... with req_pdg=%d" % int(req_pdg)

    # set ProtoShower Algo to go from data-products to a ProtoShower object
    protoshoweralg = protoshower.ProtoShowerAlgSecondShower()
    protoshoweralg.SetVertexProducer("dlsecond")
    shower_ana_unit.GetProtoShowerHelper().setProtoShowerAlg( protoshoweralg )
    
    shower_ana_unit.SetInputProducer("dlsecond")
    shower_ana_unit.SetOutputProducer("showerreco_second")

    return shower_ana_unit


