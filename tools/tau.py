import math
import ROOT

"""
Makeshift Tau object class to test out collinear mass.
"""

class Tau:

    def __init__(self,fourvect,numTrack):
        self.four_vector=fourvect
        self.pt = fourvect.Pt()
        self.px = fourvect.Px()
        self.py = fourvect.Py()
        self.eta = fourvect.Eta()
        self.phi = fourvect.Phi()
        self.E = fourvect.E()
        self.numTrack = numTrack
        self.tlv = fourvect
        self.tv2 = ROOT.TVector2()
        self.tv2.SetMagPhi(self.pt,self.phi)
