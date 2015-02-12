from ROOT import *
from tools.mass import collinearmass

import os
import sys
import math
import random
import glob

ifilepath='/cluster/data11/endw/ntuples/prod_v29/hhskim/hhskim.root'
ofilepath='/cluster/data03/sbahrase/BrtStudies/PracticeDesk/TRUTH_LEVEL_BRT_TRAINING/Full_Simulations_Samples/ggH125_mca12.root'
iFileName=ifilepath
oFileName=ofilepath

#_TreeNameInput='data12_JetTauEtmiss'

_TreeNameInput ='PowPyth8_AU2CT10_ggH125_tautauhh_mc12a'
_TreeNameOutput='Tree'
    
print "<--  input file: "+iFileName
print "--> output file: "+oFileName
iFile = TFile.Open(iFileName)
iTree = iFile.Get(_TreeNameInput)
oDirName = oFileName[::-1].split("/",1)[-1][::-1]
os.system("mkdir -p "+oDirName)
oFile = TFile(oFileName,"RECREATE")
oTree = iTree.CloneTree(0)
oTree.SetName(_TreeNameOutput)
nEntries = iTree.GetEntries()
print "number of entries in original tree = %s" %nEntries
for ientry in xrange(nEntries):
    if ientry%10000==0:
        print ientry
    ## Get the next tree in the chain and verify.                                                                                               
    if iTree.LoadTree(ientry) <  0: break
    ## Copy next entry into memory and verify.                                                                                                  
    if iTree.GetEntry(ientry) <= 0: continue
    oTree.Fill()
oTree.AutoSave()
iFile.Close()
oFile.Close()
