from larlitecv import larlitecv
from ROOT import llcv

class MatchFactory:

    def __init__(self, TYPE):
        self.types = ["shower", "track"]

        if TYPE not in self.types: 
            raise Exception("must provide : %s" % str(self.types))

        self.algo = None
        self.out_name = ""

        if TYPE == self.types[0]:
            self.algo = llcv.ShowerTruthMatch()
            self.out_name = "shower_truth_match_"

        if TYPE == self.types[1]:
            self.algo = llcv.TrackTruthMatch()
            self.out_name = "track_truth_match_"
        
        if self.algo is None:
            raise Exception("algo is none")

        if self.out_name == "":
            raise Exception("out name is empty")
        
    def get(self):
        return self.algo

    def output_file(self,num):
        return str(self.out_name + "%d.root" % num)

