# python imports
import os
import sys
import math
import random
import glob
from array import array
# ROOT imports
from ROOT import *
# local imports
from tools.mass import collinearmass
from tools import tau, BRToutput_deprecated

masses=[100,105,110,115,120,125,130,135,140,145,150]
#ifilepath='/cluster/data04/mquennev/higgs/trees/with_full_sim_variables/H/'
##ofilepath='/cluster/data04/mquennev/higgs/outputtrees/'
ifilepath="/cluster/data03/sbahrase/BrtStudies/PracticeDesk/TRUTH_LEVEL_BRT_TRAINING/Full_Simulations_Samples/"

ofilepath="/cluster/data03/sbahrase/BrtStudies/PracticeDesk/TRUTH_LEVEL_BRT_TRAINING/Full_Simulations_Samples/with_BRT_mass/VBF/"

tlv_tau1=TLorentzVector()
tlv_tau2=TLorentzVector()
tv2_met=TVector2()
for mass in masses:
    iFileName=ifilepath+'VBF_Hmass'+str(mass)+'.root'
    oFileName=ofilepath+'VBF_Hmass'+str(mass)+'.root'

    _TreeNameInput='Tree'
    _TreeNameOutput='TreeTest'
    
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

#Variables to add:
    br_alpha_vis_BRT=array('d',[0])
    br_mass_BRT=array('d',[0])

    oTree.Branch( 'alpha_vis_BRT', br_alpha_vis_BRT, 'alpha_vis_BRT/D')
    oTree.Branch( 'mass_BRT', br_mass_BRT, 'mass_BRT/D')

    for ientry in xrange(nEntries):
    ## Get the next tree in the chain and verify.
        if iTree.LoadTree(ientry) <  0: break
    ## Copy next entry into memory and verify.
        if iTree.GetEntry(ientry) <= 0: continue
        tv2_met.SetMagPhi(iTree.MET_et,iTree.MET_phi)
        tlv_tau1.SetPtEtaPhiM(iTree.tau1_pt,iTree.tau1_eta,iTree.tau1_phi,800.)
        tlv_tau2.SetPtEtaPhiM(iTree.tau2_pt,iTree.tau2_eta,iTree.tau2_phi,800.)
        tau1=tau.Tau(tlv_tau1,1)
        tau2=tau.Tau(tlv_tau2,1)
        br_alpha_vis_BRT[0]=BRToutput_deprecated.mass_BRT(tau1,tau2,tv2_met.Px(),tv2_met.Py())
        br_mass_BRT[0]=br_alpha_vis_BRT[0]*(tlv_tau1+tlv_tau2).M()
        oTree.Fill()
    oTree.AutoSave()
    iFile.Close()
    oFile.Close()
