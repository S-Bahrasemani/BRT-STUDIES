import os
from ROOT import *

print os.path.dirname(os.path.realpath(__file__))
ROOT.gROOT.LoadMacro(os.path.join(os.path.dirname(
            os.path.realpath(__file__)), "AtlasStyle.C+g"))
SetAtlasStyle()
