import ROOT
from . import tau, cutflow, variables
import math

# the place this module were called was already commented!
##import interval

def getVariableHists(treename,filelist,variabledict,recoVariableNames=True,nPoints=50,cutString='boosted'):
    
    if recoVariableNames:
        tag='_reco'
    else:
        tag='_truth'
    tlv_tau1=ROOT.TLorentzVector()
    tlv_tau2=ROOT.TLorentzVector()
    tv2_met=ROOT.TVector2()
    Tree=ROOT.TChain(treename)
    print 'Evaluating for tree, '+str(treename)+', in files: '
    for i in filelist:
        Tree.Add(i)
        print i
    nEntries=Tree.GetEntries()
    print "Entries: ", nEntries
    Hists={k:ROOT.TH1F(k+filelist[0].split('/')[-1].split('.')[0]+tag,k+filelist[0].split('/')[-1].split('.')[0]+tag,nPoints,v[-3],v[-2]) for k,v in variabledict.iteritems()}
    for ientry in xrange(nEntries):
    ## Get the next tree in the chain and verify.
        if Tree.LoadTree(ientry) <  0: break
    ## Copy next entry into memory and verify.
        if Tree.GetEntry(ientry) <= 0: continue
        if recoVariableNames:
            tv2_met.SetMagPhi(Tree.MET_et,Tree.MET_phi)
            tlv_tau1.SetPtEtaPhiM(Tree.tau1_pt,Tree.tau1_eta,Tree.tau1_phi,800.)
            tlv_tau2.SetPtEtaPhiM(Tree.tau2_pt,Tree.tau2_eta,Tree.tau2_phi,800.)
        else:
            tv2_met.SetMagPhi(Tree.met_et,Tree.met_phi)
            tlv_tau1.SetPtEtaPhiM(Tree.lep1_pt,Tree.lep1_eta,Tree.lep1_phi,800.)
            tlv_tau2.SetPtEtaPhiM(Tree.lep2_pt,Tree.lep2_eta,Tree.lep2_phi,800.)
        tau1=tau.Tau(tlv_tau1,1)
        tau2=tau.Tau(tlv_tau2,1)
        if recoVariableNames: 
            if not(cutflow.passedCut(Tree,cutString)):
                continue
        else:
            if not(cutflow.passedCut_truth(Tree)):
                continue
        for k,v in variabledict.iteritems():
            if not(k=='met_et'):
                Hists[k].Fill(variables.evalVariable(k,tau1,tau2,tv2_met.Px(),tv2_met.Py()))
            else:
                Hists[k].Fill(variables.evalVariable(k,tau1,tau2,tv2_met.Px(),tv2_met.Py())/1000.)
    return Hists


def doZCal(treename,filelist,reader,nPoints=100,cutString='preselection',Zmass=91.1876):
    tlv_tau1=ROOT.TLorentzVector()
    tlv_tau2=ROOT.TLorentzVector()
    tv2_met=ROOT.TVector2()
    Tree=ROOT.TChain(treename)
    print 'Evaluating for tree, '+str(treename)+', in files: '
    for i in filelist:
        Tree.Add(i)
        print i
    nEntries=Tree.GetEntries()
    print "Entries: ", nEntries
    mass_BRT=ROOT.TH1F('massBRT'+filelist[0].split('/')[-1].split('.')[0],'massBRT',nPoints,0,250.)

    for ientry in xrange(nEntries):
    ## Get the next tree in the chain and verify.
        if Tree.LoadTree(ientry) <  0: break
    ## Copy next entry into memory and verify.
        if Tree.GetEntry(ientry) <= 0: continue
        tv2_met.SetMagPhi(Tree.met_et,Tree.met_phi)
        tlv_tau1.SetPtEtaPhiM(Tree.lep1_pt,Tree.lep1_eta,Tree.lep1_phi,800.)
        tlv_tau2.SetPtEtaPhiM(Tree.lep2_pt,Tree.lep2_eta,Tree.lep2_phi,800.)
        tau1=tau.Tau(tlv_tau1,1)
        tau2=tau.Tau(tlv_tau2,1)
        if not(cutflow.passedCut_truth(Tree)):
            continue
        brt=reader.mass_BRT(tau1,tau2,tv2_met.Px(),tv2_met.Py())
        mass_BRT.Fill(brt)
    #hist_interval=interval.getInterval(mass_BRT,0.8)
    #mass_BRT.Fit('gaus','rs','',hist_interval[0],hist_interval[1])
    #fit_func=mass_BRT.GetFunction('gaus')
    #mean=fit_func.GetParameter(1)
    return Zmass-mass_BRT.GetMean(),mass_BRT

def getHists(treename,filename,reader,lumi,overallLumi,nPoints=100,cutString='boosted',ZCalConst=0.,tauVeto=False):
    tlv_tau1=ROOT.TLorentzVector()
    tlv_tau2=ROOT.TLorentzVector()
    tv2_met=ROOT.TVector2()
    Tree=ROOT.TChain(treename)
    print 'Evaluating for tree, '+str(treename)+', in file: '+filename
    Tree.Add(filename)
    nEntries=Tree.GetEntries()
    print "Entries: ", nEntries
    mass_BRT=ROOT.TH1F('massBRT'+filename,'massBRT',nPoints,0,250.)
    mass_MMC=ROOT.TH1F('massMMC'+filename,'massMMC',nPoints,0,250.)
    mass_Col=ROOT.TH1F('massCol'+filename,'massCol',nPoints,0,250.)
    mass_BRT.Sumw2()
    mass_MMC.Sumw2()
    mass_Col.Sumw2()
    for ientry in xrange(nEntries):
    ## Get the next tree in the chain and verify.
        if Tree.LoadTree(ientry) <  0: break
    ## Copy next entry into memory and verify.
        if Tree.GetEntry(ientry) <= 0: continue
        tv2_met.SetMagPhi(Tree.MET_et,Tree.MET_phi)
        tlv_tau1.SetPtEtaPhiM(Tree.tau1_pt,Tree.tau1_eta,Tree.tau1_phi,800.)
        tlv_tau2.SetPtEtaPhiM(Tree.tau2_pt,Tree.tau2_eta,Tree.tau2_phi,800.)
        tau1=tau.Tau(tlv_tau1,1)
        tau2=tau.Tau(tlv_tau2,1)
        if not(cutflow.passedCut(Tree,cutString)):
            continue
        if tauVeto and event.tau1_charge*event.tau2_charge==-1:
            continue
        brt=reader.mass_BRT(tau1,tau2,tv2_met.Px(),tv2_met.Py())+ZCalConst
        mass_BRT.Fill(brt,overallLumi/lumi)
        
        if Tree.mmc1_resonance_m>0.:
            mass_MMC.Fill(Tree.mmc1_resonance_m,overallLumi/lumi)
        if  Tree.mass_collinear_tau1_tau2>0.:
            mass_Col.Fill(Tree.mass_collinear_tau1_tau2/1000.,overallLumi/lumi)
    return [mass_BRT,mass_MMC,mass_Col]

def getPileupDependence(treename,filename,reader,truthMass,lumi,OverallLumi,nPoints=100,cutString='boosted',ZCalConst=0.):
    tlv_tau1=ROOT.TLorentzVector()
    tlv_tau2=ROOT.TLorentzVector()
    tv2_met=ROOT.TVector2()
    Tree=ROOT.TChain(treename)
    print 'Evaluating for tree, '+str(treename)+', in file: '+filename
    Tree.Add(filename)
    nEntries=Tree.GetEntries()
    print "Entries: ", nEntries
    mass_BRT={}
    mass_MMC={}
    mass_Col={}

    mass_BRT_hist=ROOT.TH1F('massBRT'+filelist[0].split('/')[-1].split('.')[0],'massBRT',nPoints,0,250.)
    mass_MMC_hist=ROOT.TH1F('massMMC'+filelist[0].split('/')[-1].split('.')[0],'massMMC',nPoints,0,250.)
    mass_Col_hist=ROOT.TH1F('massCol'+filelist[0].split('/')[-1].split('.')[0],'massCol',nPoints,0,250.)
    for ientry in xrange(nEntries):
    ## Get the next tree in the chain and verify.
        if Tree.LoadTree(ientry) <  0: break
    ## Copy next entry into memory and verify.
        if Tree.GetEntry(ientry) <= 0: continue
        tv2_met.SetMagPhi(Tree.MET_et,Tree.MET_phi)
        tlv_tau1.SetPtEtaPhiM(Tree.tau1_pt,Tree.tau1_eta,Tree.tau1_phi,800.)
        tlv_tau2.SetPtEtaPhiM(Tree.tau2_pt,Tree.tau2_eta,Tree.tau2_phi,800.)
        pileup=Tree.averageIntPerXing
        tau1=tau.Tau(tlv_tau1,1)
        tau2=tau.Tau(tlv_tau2,1)
        if not(cutflow.passedCut(Tree,cutString)):
            continue
        brt=reader.mass_BRT(tau1,tau2,tv2_met.Px(),tv2_met.Py())+ZCalConst
        if pileup in mass_BRT.keys():
            mass_BRT[pileup].Fill(brt)
            mass_MMC[pileup].Fill(Tree.mmc1_resonance_m)
            mass_Col[pileup].Fill(Tree.mass_collinear_tau1_tau2/1000.)
        else:
            mass_BRT[pileup]=ROOT.TH1F('massBRT'+filelist[0].split('/')[-1].split('.')[0]+'_'+str(pileup),'massBRT',nPoints,0,250.)
            mass_MMC[pileup]=ROOT.TH1F('massMMC'+filelist[0].split('/')[-1].split('.')[0]+'_'+str(pileup),'massMMC',nPoints,0,250.)
            mass_Col[pileup]=ROOT.TH1F('massCol'+filelist[0].split('/')[-1].split('.')[0]+'_'+str(pileup),'massCol',nPoints,0,250.)
            
            mass_BRT[pileup].Fill(brt)
            mass_MMC[pileup].Fill(Tree.mmc1_resonance_m)
            mass_Col[pileup].Fill(Tree.mass_collinear_tau1_tau2/1000.)
        mass_BRT_hist.Fill(brt)
        mass_MMC_hist.Fill(Tree.mmc1_resonance_m)
        mass_Col_hist.Fill(Tree.mass_collinear_tau1_tau2/1000.)

        graphs_BRT=ROOT.TGraphErrors(len(mass_BRT.keys()))
        graphs_MMC=ROOT.TGraphErrors(len(mass_MMC.keys()))
        graphs_Col=ROOT.TGraphErrors(len(mass_Col.keys()))
        for ikey,key in enumerate(mass_BRT.keys()):
            graphs_BRT.SetPoint(ikey,float(key),mass_BRT[key].GetMean()-truthMass)
            graphs_MMC.SetPoint(ikey,float(key),mass_MMC[key].GetMean()-truthMass)
            graphs_Col.SetPoint(ikey,float(key),mass_Col[key].GetMean()-truthMass)

            graphs_BRT.SetPointError(ikey,0,mass_BRT[key].GetRMS()/math.sqrt(mass_BRT[key].GetEntries()))
            graphs_MMC.SetPointError(ikey,0,mass_MMC[key].GetRMS()/math.sqrt(mass_MMC[key].GetEntries()))
            graphs_Col.SetPointError(ikey,0,mass_Col[key].GetRMS()/math.sqrt(mass_Col[key].GetEntries()))
    mass_BRT_hist.Sumw2()
    mass_MMC_hist.Sumw2()
    mass_Col_hist.Sumw2()
    mass_BRT_hist.Scale(overallLumi/lumi)
    mass_MMC_hist.Scale(overallLumi/lumi)
    mass_Col_hist.Scale(overallLumi/lumi)
    return [graphs_BRT,graphs_MMC,graphs_Col],[mass_BRT_hist,mass_MMC_hist,mass_Col_hist]
