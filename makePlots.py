import math
import array
# import interval

import ROOT

from tools import testOnTree, BRToutput, roccurve
import tools.atlasstyle

############################################################################
##############  Configure the options in this section  #####################
#### Currently, variable list must be updated in ./modules/BRTOutput.py ####

### Choose analysis to run (One or both)###
_compareMassEstimators=True
_variablePlots=False

### Use which samples (_compareMassEstimators only) ###
_useH=True
_useZ=True
_doZCal=False

### Fraction of integral of histogram to do gaussian fits over ###
_integralFraction=0.8

### Where to send plots ###
_outputFilePath='./plots'

### Full sim sample settings (_compareMassEstimators only) ###
_overallLumi=20.3

#Z+jets samples
_pathToZ='/cluster/data03/sbahrase/BrtStudies/PracticeDesk/TRUTH_LEVEL_BRT_TRAINING/Full_Simulations_Samples/'
_fileListZ=[_pathToZ+'ZtautauNp'+str(i)+'.root' for i in range(6)]
_treeNameZ='Tree'
_ZLumiList=[22.8,51.4,54.3,59.7,84.4,155.1]

#DY Z+jets samples
_pathToDYZ='/cluster/data03/sbahrase/BrtStudies/PracticeDesk/TRUTH_LEVEL_BRT_TRAINING/Full_Simulations_Samples/'
_fileListDYZ=[_pathToDYZ+'DY_ZtautauNp'+str(i)+'.root' for i in range(6)]
_treeNameDYZ='Tree'
_DYZLumiList=[.2,2.3,7.4,10.7,119.9,145.5]

#EW Z+jets samples
_pathToEWZ='/cluster/data03/sbahrase/BrtStudies/PracticeDesk/TRUTH_LEVEL_BRT_TRAINING/Full_Simulations_Samples/'
_fileListEWZ=[_pathToEWZ+'EW_Ztautau.root' for i in range(1)]
_treeNameEWZ='Tree'
_EWZLumiList=[1390.0]

#Higgs
_masses=[100,105,110,115,120,125,130,135,140,145,150]

#ggH
_pathToggH='/cluster/data03/sbahrase/BrtStudies/PracticeDesk/TRUTH_LEVEL_BRT_TRAINING/Full_Simulations_Samples/'
_fileListggH=[_pathToggH+'ggH_Hmass'+str(i)+'.root' for i in _masses]
_treeNameggH='Tree'
_ggHLumiList=[1454.0,1619.3,1819.7,2080.3,2433.2,3912.1,3673.1,4802.9,6144.4,9468.1,14716.7]

#VBFH
_pathToVBFH='/cluster/data03/sbahrase/BrtStudies/PracticeDesk/TRUTH_LEVEL_BRT_TRAINING/Full_Simulations_Samples/'
_fileListVBFH=[_pathToVBFH+'VBFH_Hmass'+str(i)+'.root' for i in _masses]
_treeNameVBFH='Tree'
_VBFHLumiList=[18813.9,19981.9,21451.0,23629.9,26676.8,47773.8,37606.8,47633.5,58468.3,88671.4,135924.5]

### Truth Level samples (for _doZCal or _variablePlots) ###
#Truth Higgs
_pathToTruth='/media/Portable Drive/Work/higgs/rootfiles/'
_fileListTruth=[_pathToTruth+'Hmass'+str(i)+'.root' for i in _masses]
_treeNameTruth='Tree'
#Truth Z
_fileListTruthZ=[_pathToTruth+'Z.root']
_treeNameTruthZ='Tree'

### Apply which cuts? (preselection, boosted, vbf) ###
_cutString='boosted'

### Are energy variables scaled to be relative? To which variable? (If not sure, use False) ###
_useRelVariables=False
_EnergyScale='ptsum_lep1_lep2_met'

### Is target variable a correction factor? To which variable? (If not sure use, False) ###
_isCorrection=False
_CorrectionFactor='visible_mass'

### BRT weight file ###
##_weightFilePath='/media/Portable Drive/Work/higgs/weights/applying_preselection_cuts/with_analysis_variables/hadhad/TMVARegression_BRT_HiggsMass_nEventsMin40_AdaBoostBeta0.2_nCuts20_MaxDepth100_NTrees80.weights.xml'

_weightFilePath='/cluster/data03/sbahrase/BrtStudies/PracticeDesk/TRUTH_LEVEL_BRT_TRAINING/weights/without_preselection_cuts/with_analysis_variables/hadhad/TMVARegression_BRT_HiggsMass_nEventsMin20_AdaBoostBeta0.2_nCuts20_MaxDepth100_NTrees40.weights.xml'

 
#_weightFilePath='/media/Portable Drive/Work/higgs/weights/applying_preselection_cuts/with_analysis_variables/hadhad/TMVARegression_BRT_HiggsMass_nEventsMin40_AdaBoostBeta0.2_nCuts20_MaxDepth100_NTrees80.weights.xml'
### Number of bins in mass distributions for calculations ###
_nPoints=100

### Number of bins in mass plots (Fewer than or equal to _nPoints) ###
_nBins=50

### Variable dictionary. (_variablePlots only) ###
variabledict={}
variabledict['lep1_pt']                      = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,200000., 40]
variabledict['lep1_eta']                     = [''   ,'F',-10.00,10.00,  array.array('f',[0]),ROOT.TBranch(), -2.5,2.5,  200]
variabledict['lep2_pt']                      = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,200000., 40]
variabledict['lep2_eta']                     = [''   ,'F',-10.00,10.00,  array.array('f',[0]),ROOT.TBranch(), -2.5,2.5,  200]
variabledict['met_et']                       = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,200., 20]
variabledict['transverse_mass_lep1_lep2']    = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
variabledict['transverse_mass_lep1_met']     = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
variabledict['transverse_mass_lep2_met']     = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
variabledict['dphi_lep1_met']                = ['rad','F',  0.00, 3.15,  array.array('f',[0]),ROOT.TBranch(),   0.00, 3.15,   30]
variabledict['dphi_lep2_met']                = ['rad','F',  0.00, 3.15,  array.array('f',[0]),ROOT.TBranch(),   0.00, 3.15,   30]
variabledict['dphi_lep_lep']                 = ['rad','F',  0.00, 3.15,  array.array('f',[0]),ROOT.TBranch(),   0.00, 3.15,   30]
variabledict['deta_lep_lep']                 = ['',   'F',  0.00,20.00,  array.array('f',[0]),ROOT.TBranch(),   0.00,3.15,  200]
variabledict['dR_lep_lep']                   = ['',   'F',  0.00,25.00,  array.array('f',[0]),ROOT.TBranch(),   0.00,3.15,  200]
variabledict['ptsum_lep1_lep2_met']          = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,500000., 30]
variabledict['ptsum_lep1_lep2']              = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,300000., 25]
variabledict['pttot_lep1_lep2_met']          = ['',   'F',  0.00, 2.00,  array.array('f',[0]),ROOT.TBranch(),   0.00, 1.10,   22]
variabledict['pttot_lep1_lep2']              = ['',   'F',  0.00, 2.00,  array.array('f',[0]),ROOT.TBranch(),   0.00, 1.10,   22]
variabledict['ptdiff_lep1_lep2']             = ['',   'F',  0.00, 2.00,  array.array('f',[0]),ROOT.TBranch(),   0.00, 1.10,   22]
variabledict['visible_mass']                 = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
variabledict['collinear_mass']               = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
variabledict['proj_mass_lep1']               = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
variabledict['proj_mass_lep2']               = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
variabledict['met_phi_centrality']           = ['',   'F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   -1.45,  1.45, 50]

############################## Should be no need to edit below here ##############################
##################################################################################################

with open(_outputFilePath+'/output.txt', 'w') as f:
    f.write('Analysis Output.\n')
    f.write('Weight File:\t'+_weightFilePath+'\n')
    
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetFillColor(0)

_Zmass=91.1876
if _compareMassEstimators:
    reader=BRToutput.BRTMass(_weightFilePath)

    if _useRelVariables:
        reader.useRelativeVariables(_EnergyScale)

    if _isCorrection:
        reader.useCorrectionAsTarget(_CorrectionFactor)
else:
    _useH=False
    _useZ=False
    _doZCal=False

if _variablePlots:
    print 'Calculating variables in Truth Trees...'
    variable_hists_Truth={imass:testOnTree.getVariableHists(_treeNameTruth,[i],variabledict,recoVariableNames=False) for imass,i in zip(_masses,_fileListTruth)}
    print 'Complete!\n'
    print 'Calculating variables in Reco Trees...'
    variable_hists_H={imass:testOnTree.getVariableHists(_treeNameggH,[i],variabledict,recoVariableNames=True) for imass,i in zip(_masses,_fileListggH)}
    print 'Complete!\n'

if _doZCal:
    print "Calibrating BRT on Z tree..."
    const,ztruthhist=testOnTree.doZCal(_treeNameTruthZ,_fileListTruthZ,reader,nPoints=_nPoints,Zmass=_Zmass)
    print "Offset = ", const
    print "Complete!\n"
else:
    const=0.

if _useH:
    print 'Evaluating on H samples...'
    H_hists={imass:testOnTree.getHists(_treeNameggH,i,reader,ilumi,_overallLumi,nPoints=_nPoints,cutString=_cutString,ZCalConst=const) for imass,i,ilumi in zip(_masses,_fileListggH,_ggHLumiList)}
    H_VBFHhists={imass:testOnTree.getHists(_treeNameVBFH,i,reader,ilumi,_overallLumi,nPoints=_nPoints,cutString=_cutString,ZCalConst=const) for imass,i,ilumi in zip(_masses,_fileListVBFH,_VBFHLumiList)}
    for imass in H_hists.keys():
        for i in range(3):
            H_hists[imass][i].Add(H_VBFHhists[imass][i])
    print 'Complete!\n'

if _useZ:
    print 'Evaluating on Z samples...'
    Z_hists=testOnTree.getHists(_treeNameEWZ,_fileListEWZ[0],reader,_ZLumiList[0],_overallLumi,nPoints=_nPoints,cutString=_cutString,ZCalConst=const)
    for np in range(6):
        newZhists=testOnTree.getHists(_treeNameZ,_fileListZ[np],reader,_ZLumiList[np],_overallLumi,nPoints=_nPoints,cutString=_cutString,ZCalConst=const)
        newDYZhists=testOnTree.getHists(_treeNameDYZ,_fileListDYZ[np],reader,_DYZLumiList[np],_overallLumi,nPoints=_nPoints,cutString=_cutString,ZCalConst=const)
        for i in range(3):
            Z_hists[i].Add(newZhists[i])
            Z_hists[i].Add(newDYZhists[i])
    print 'Complete!\n'

print 'Making Plots...'
canvas=ROOT.TCanvas()

dist_label_1=ROOT.TLatex(0.2,0.85,"#splitline{H(125)#rightarrow#tau_{had}#tau_{had}}{Boosted}")
dist_label_1.SetNDC()
dist_label_2=ROOT.TLatex(0.65,0.6,"#splitline{#it{#bf{ATLAS}}}{Work in Progress}")
dist_label_2.SetNDC()
reco_vs_truth_label=ROOT.TLatex(0.55,0.25,"#splitline{H(125)#rightarrow#tau_{had}#tau_{had} Boosted}{#it{#bf{ATLAS}} Work in Progress}")
reco_vs_truth_label.SetNDC()
rms_label=ROOT.TLatex(0.2,0.8,"#splitline{H(125)#rightarrow#tau_{had}#tau_{had} Boosted}{#it{#bf{ATLAS}} Work in Progress}")
rms_label.SetNDC()
roc_label=ROOT.TLatex(0.2,0.45,"#splitline{H(125)#rightarrow#tau_{had}#tau_{had} Boosted}{#it{#bf{ATLAS}} Work in Progress}")
roc_label.SetNDC()

#========== Z mass distribution ==========#
if _useZ:
    Zstack=ROOT.THStack()
    Zlegend=ROOT.TLegend(0.7,0.7,0.9,0.9)

    Z_hists[0].SetLineColor(ROOT.kRed)
    Z_hists[0].SetLineWidth(2)
    Zstack.Add(Z_hists[0].Rebin(_nPoints/_nBins,"brt_z_hist"))
    Zlegend.AddEntry(Z_hists[0],'BRT mass','l')

    Z_hists[1].SetLineColor(ROOT.kBlack)
    Z_hists[1].SetLineWidth(2)
    Zstack.Add(Z_hists[1].Rebin(_nPoints/_nBins,"mmc_z_hist"))
    Zlegend.AddEntry(Z_hists[1],'MMC mass','l')
    
    Z_hists[2].SetLineColor(ROOT.kBlue)
    Z_hists[2].SetLineWidth(2)
    Zstack.Add(Z_hists[2].Rebin(_nPoints/_nBins,"col_z_hist"))
    Zlegend.AddEntry(Z_hists[2],'Collinear mass','l')
    
    Zstack.Draw('NOSTACK hist')

    Zstack.SetTitle('Reconstructed Z Mass Distribution')
    Zstack.GetXaxis().SetTitle('Reconstructed Mass (GeV)')
    Zstack.GetYaxis().SetTitle('Events / 5 GeV')

    Zlegend.Draw('same')
    Zlegend.SetTextSize(0.0375)
    dist_label_1=ROOT.TLatex(0.2,0.85,"#splitline{Z#rightarrow#tau_{had}#tau_{had}}{Boosted}")
    dist_label_1.SetNDC()
    dist_label_2=ROOT.TLatex(0.65,0.6,"#splitline{#it{#bf{ATLAS}}}{Work in Progress}")
    dist_label_2.SetNDC()
    dist_label_1.Draw('same')
    dist_label_2.Draw('same')
    canvas.SaveAs(_outputFilePath+'/Z_mass.png')
    canvas.Clear()

#========== Truth Z mass distribution ==========#
if _doZCal:
    ztruthhist.SetLineColor(ROOT.kRed)
    ztruthhist.SetLineWidth(2)
    ztruthhist.Rebin(_nPoints/_nBins)
    ztruthhist.Scale(1/ztruthhist.GetEntries())
    Ztruthlegend=ROOT.TLegend(0.7,0.7,0.9,0.9)
    ztruthhist.Draw('hist')
    truemass=ROOT.TLine(91.1876,0,91.1876,ztruthhist.GetMaximum())
    truemass.SetLineWidth(2)
    truemass.SetLineStyle(ROOT.kDashed)
    truemass.Draw('same')
    Ztruthlegend.AddEntry(ztruthhist,'BRT Mass','l')
    Ztruthlegend.AddEntry(truemass,'True Z Mass','l')
    ztruthhist.GetXaxis().SetTitle('Reconstructed Mass (GeV)')
    ztruthhist.GetYaxis().SetTitle('Events / 5 GeV')
    Ztruthlegend.Draw('same')
    Ztruthlegend.SetTextSize(0.0375)
    dist_label_1=ROOT.TLatex(0.2,0.85,"#splitline{Truth Z#rightarrow#tau_{had}#tau_{had}}{Preselection}")
    dist_label_1.SetNDC()
    dist_label_2=ROOT.TLatex(0.65,0.6,"#splitline{#it{#bf{ATLAS}}}{Work in Progress}")
    dist_label_2.SetNDC()
    dist_label_1.Draw('same')
    dist_label_2.Draw('same')
    canvas.SaveAs(_outputFilePath+'/Z_mass_truth.png')
    canvas.Clear()

#========== H mass distributions ==========#
if _useH:
    for i in _masses:
        Hstack=ROOT.THStack()
        Hlegend=ROOT.TLegend(0.7,0.7,0.9,0.9)

        H_hists[i][0].SetLineColor(ROOT.kRed)
        H_hists[i][0].SetLineWidth(2)
        Hstack.Add(H_hists[i][0].Rebin(_nPoints/_nBins,"brt_h_"+str(i)+"_hist"))
        Hlegend.AddEntry(H_hists[i][0],'BRT mass','l')

        H_hists[i][1].SetLineColor(ROOT.kBlack)
        H_hists[i][1].SetLineWidth(2)
        Hstack.Add(H_hists[i][1].Rebin(_nPoints/_nBins,"mmc_h_"+str(i)+"_hist"))
        Hlegend.AddEntry(H_hists[i][1],'MMC mass','l')
        
        H_hists[i][2].SetLineColor(ROOT.kBlue)
        H_hists[i][2].SetLineWidth(2)
        Hstack.Add(H_hists[i][2].Rebin(_nPoints/_nBins,"col_h_"+str(i)+"_hist"))
        Hlegend.AddEntry(H_hists[i][2],'Collinear mass','l')
        
        Hstack.Draw('NOSTACK hist')

        Hstack.SetTitle('Reconstructed H('+str(i)+') Mass Distribution')
        Hstack.GetXaxis().SetTitle('Reconstructed Mass (GeV)')
        Hstack.GetYaxis().SetTitle('Events / 5 GeV')

        Hlegend.Draw('same')
        Hlegend.SetTextSize(0.0375)
        dist_label_1=ROOT.TLatex(0.2,0.85,"#splitline{H("+str(i)+")#rightarrow#tau_{had}#tau_{had}}{Boosted}")
        dist_label_1.SetNDC()
        if i<=125:
            dist_label_2=ROOT.TLatex(0.65,0.6,"#splitline{#it{#bf{ATLAS}}}{Work in Progress}")
        else:
            dist_label_2=ROOT.TLatex(0.2,0.6,"#splitline{#it{#bf{ATLAS}}}{Work in Progress}")
        dist_label_2.SetNDC()
        dist_label_1.Draw('same')
        dist_label_2.Draw('same')
        canvas.SaveAs(_outputFilePath+'/H_mass_'+str(i)+'.png')
        canvas.Clear()

#========== H_125 vs Z mass distribution ==========#
if _useH and _useZ:
    Z_hists_norm=[i.Clone() for i in Z_hists]
    Z_hists_norm=[i.Scale(1./i.GetEntries()) for i in Z_hists_norm]

    H_hists_norm={k:[i.Clone() for i in v] for k,v in H_hists.iteritems()}
    H_hists_norm={k:[i.Scale(1./i.GetEntries()) for i in v] for k,v in H_hists_norm.iteritems()}

    ZHstack=[ROOT.THStack() for i in range(3)]

    for i in range(3):
        if i==0:
            tag='BRT'
        elif i==1:
            tag='MMC'
        else:
            tag='Collinear'

        ZHstack[i].Add(Z_hists[i].Rebin(_nPoints/_nBins,"z_"+str(i)+"_hist"))
        ZHstack[i].Add(H_hists[125][i].Rebin(_nPoints/_nBins,"h_125"+str(i)+"_hist"))

        ZHstack[i].Draw('nostack hist')
        dist_label_1=ROOT.TLatex(0.2,0.85,"#splitline{Z/H#rightarrow#tau_{had}#tau_{had}}{Boosted}")
        dist_label_1.SetNDC()
        dist_label_2=ROOT.TLatex(0.65,0.6,"#splitline{#it{#bf{ATLAS}}}{Work in Progress}")
        dist_label_2.SetNDC()
        dist_label_1.Draw('same')
        dist_label_2.Draw('same')
        ZHstack[i].SetTitle('Reconstructed Z and H(125) Mass Distribution ('+tag+')')
        ZHstack[i].GetXaxis().SetTitle('Reconstructed Mass (GeV)')
        ZHstack[i].GetYaxis().SetTitle('Entries')

        if i==0:
            canvas.SaveAs(_outputFilePath+'/Z_H_mass_brt.png')
        elif i==1:
            canvas.SaveAs(_outputFilePath+'/Z_H_mass_mmc.png')
        else:
            canvas.SaveAs(_outputFilePath+'/Z_H_mass_col.png')

        canvas.Clear()

#========== Output vs. Truth mass ==========#
if _useH:
    nPoints=len(_masses)

    if _useZ:
        nPoints+=1

    meanlegend=ROOT.TLegend(0.2,0.7,0.4,0.9)
    reslegend=ROOT.TLegend(0.2,0.7,0.4,0.9)
    
    meanGraphs=[ROOT.TGraphErrors(nPoints) for i in range(3)]
    rmsGraphs=[ROOT.TGraph(nPoints) for i in range(3)]
    resGraphs=[ROOT.TGraph(nPoints) for i in range(3)]

    masslist=_masses
    if _useZ:
        masslist=[_Zmass]+masslist
    Z_mean=[]
    Z_rms=[]
    H_mean={}
    H_rms={}
    for i,imass in enumerate(masslist):
        H_mean[imass]=[]
        H_rms[imass]=[]

        for j in range(3):
            if j==0:
                offset=-1.
            elif j==1:
                offset=1.
            else:
                offset=0

            if i==0 and _useZ:
                #hist_interval=interval.getInterval(Z_hists[j],_integralFraction)
                #Z_hists[j].Fit('gaus','rs','',hist_interval[0],hist_interval[1])
                #fit_func=Z_hists[j].GetFunction('gaus')
                Z_mean.append(Z_hists[j].GetMean())#Z_hists[j].GetMean())#fit_func.GetParameter(1))
                Z_rms.append(Z_hists[j].GetRMS())#Z_hists[j].GetRMS())#fit_func.GetParameter(2))
                meanGraphs[j].SetPoint(i,float(imass+offset),Z_mean[j])
                meanGraphs[j].SetPointError(i,0,Z_rms[j])
                rmsGraphs[j].SetPoint(i,float(imass),Z_rms[j])
                resGraphs[j].SetPoint(i,float(imass),Z_rms[j]/float(imass))
            else:
                #hist_interval=interval.getInterval(H_hists[imass][j],_integralFraction)
                #H_hists[imass][j].Fit('gaus','rs','',hist_interval[0],hist_interval[1])
                #fit_func=H_hists[imass][j].GetFunction('gaus')
                H_mean[imass].append(H_hists[imass][j].GetMean())#H_hists[imass][j].GetMean())#fit_func.GetParameter(1))
                H_rms[imass].append(H_hists[imass][j].GetRMS())#H_hists[imass][j].GetRMS())#fit_func.GetParameter(2))
                meanGraphs[j].SetPoint(i,float(imass+offset),H_mean[imass][j])
                meanGraphs[j].SetPointError(i,0,H_rms[imass][j])
                rmsGraphs[j].SetPoint(i,float(imass),H_rms[imass][j])
                resGraphs[j].SetPoint(i,float(imass),H_rms[imass][j]/float(imass))

    meanMultiGraph=ROOT.TMultiGraph()
    rmsMultiGraph=ROOT.TMultiGraph()
    resMultiGraph=ROOT.TMultiGraph()

    for i in range(3):
        if i==0:
            color=ROOT.kRed    
            title='BRT mass'
        elif i==1:
            color=ROOT.kBlack
            title='MMC mass'
        else:
            color=ROOT.kBlue
            title='Collinear mass'

        meanGraphs[i].SetMarkerColor(color)
        meanGraphs[i].SetLineColor(color)
        rmsGraphs[i].SetMarkerColor(color)
        resGraphs[i].SetMarkerColor(color)
        meanGraphs[i].SetMarkerStyle(2)
        rmsGraphs[i].SetMarkerStyle(8)
        resGraphs[i].SetMarkerStyle(8)
        meanGraphs[i].SetLineWidth(2)

        
        meanlegend.AddEntry(meanGraphs[i],title,'P')
        reslegend.AddEntry(resGraphs[i],title,'P')
        meanMultiGraph.Add(meanGraphs[i])
        rmsMultiGraph.Add(rmsGraphs[i])
        resMultiGraph.Add(resGraphs[i])

    meanMultiGraph.Draw('AP')

    meanMultiGraph.SetTitle('Reconstructed Mass vs. Truth Mass')

    line=ROOT.TLine(meanMultiGraph.GetXaxis().GetXmin(),meanMultiGraph.GetXaxis().GetXmin(),meanMultiGraph.GetXaxis().GetXmax(),meanMultiGraph.GetXaxis().GetXmax())
    line.SetLineWidth(2)

    meanlegend.AddEntry(line,"m_{Reco} = m_{True}",'L')

    line.Draw('same')

    meanlegend.Draw('same')
    meanlegend.SetTextSize(0.0375)
    reco_vs_truth_label.Draw('same')
    meanMultiGraph.SetTitle('Reconstructed Mass vs. Truth Mass')
    meanMultiGraph.GetXaxis().SetTitle('Truth Mass (GeV)')
    meanMultiGraph.GetYaxis().SetTitle('Reconstructed Mass (GeV)')

    canvas.SaveAs(_outputFilePath+'/reco_mass_vs_truth.png')
    canvas.Clear()

    rmsMultiGraph.Draw('AP')
    Hlegend.Draw('same')
    Hlegend.SetTextSize(0.0375)
    rms_label.Draw('same')
    rmsMultiGraph.SetMinimum(0.)
    rmsMultiGraph.SetMaximum(50.)

    rmsMultiGraph.SetTitle('Reconstructed Mass RMS vs. Truth Mass')
    rmsMultiGraph.GetXaxis().SetTitle('Truth Mass (GeV)')
    rmsMultiGraph.GetYaxis().SetTitle('Reco Mass RMS (GeV)')

    canvas.SaveAs(_outputFilePath+'/reco_rms_vs_truth.png')
    canvas.Clear()

    resMultiGraph.Draw('AP')

    Hlegend.Draw('same')
    Hlegend.SetTextSize(0.0375)
    rms_label.Draw('same')
    resMultiGraph.SetMinimum(0.)
    resMultiGraph.SetMaximum(0.5)

    resMultiGraph.SetTitle('Reconstructed Mass Resolution vs. Truth Mass')
    resMultiGraph.GetXaxis().SetTitle('Truth Mass (GeV)')
    resMultiGraph.GetYaxis().SetTitle('Reco Mass Resolution')

    canvas.SaveAs(_outputFilePath+'/reco_res_vs_truth.png')
    canvas.Clear()

#========== Z vs. H125 ROC ==========#
if _useH and _useZ:
    roc_curves_Z_125=[roccurve.roccurve(H_hists[125][i],Z_hists[i]) for i in range(3)]
    roc_legend=ROOT.TLegend(0.2,0.2,0.4,0.4)
    
    roc_curves_Z_125[0].SetLineWidth(2)
    roc_curves_Z_125[0].SetLineColor(ROOT.kRed)
    roc_legend.AddEntry(roc_curves_Z_125[0],'BRT Mass','L')
    
    roc_curves_Z_125[1].SetLineWidth(2)
    roc_curves_Z_125[1].SetLineColor(ROOT.kBlack)
    roc_legend.AddEntry(roc_curves_Z_125[1],'MMC Mass','L')
    
    roc_curves_Z_125[2].SetLineWidth(2)
    roc_curves_Z_125[2].SetLineColor(ROOT.kBlue)
    roc_legend.AddEntry(roc_curves_Z_125[2],'Collinear Mass','L')
    
    for i in range(3):
        if i==0:
            roc_curves_Z_125[i].Draw('AC')
            roc_curves_Z_125[i].GetXaxis().SetRangeUser(0.,1.)
            roc_curves_Z_125[i].GetYaxis().SetRangeUser(0.,1.)
            roc_curves_Z_125[i].SetTitle('Z vs H(125) ROC Curve')
            roc_curves_Z_125[i].GetXaxis().SetTitle('H(125) Acceptance')
            roc_curves_Z_125[i].GetYaxis().SetTitle('Z Rejection')
        else:
            roc_curves_Z_125[i].Draw('C')
    roc_legend.Draw('same')
    roc_legend.SetTextSize(0.0375)
    roc_label.Draw('same')
    canvas.SaveAs(_outputFilePath+'/roc.png')

#========== H120 vs. H125 ROC ==========#
if _useH:
    roc_curves_120_125=[roccurve.roccurve(H_hists[125][i],H_hists[120][i]) for i in range(3)]
    roc_legend=ROOT.TLegend(0.2,0.2,0.4,0.4)
    
    roc_curves_120_125[0].SetLineWidth(2)
    roc_curves_120_125[0].SetLineColor(ROOT.kRed)
    roc_legend.AddEntry(roc_curves_120_125[0],'BRT Mass','L')
    
    roc_curves_120_125[1].SetLineWidth(2)
    roc_curves_120_125[1].SetLineColor(ROOT.kBlack)
    roc_legend.AddEntry(roc_curves_120_125[1],'MMC Mass','L')

    roc_curves_120_125[2].SetLineWidth(2)
    roc_curves_120_125[2].SetLineColor(ROOT.kBlue)
    roc_legend.AddEntry(roc_curves_120_125[2],'Collinear Mass','L')
    
    for i in range(3):
        if i==0:
            roc_curves_120_125[i].Draw('AC')
            roc_curves_120_125[i].GetXaxis().SetRangeUser(0.,1.)
            roc_curves_120_125[i].GetYaxis().SetRangeUser(0.,1.)
            roc_curves_120_125[i].SetTitle('H(120) vs H(125) ROC Curve')
            roc_curves_120_125[i].GetXaxis().SetTitle('H(125) Acceptance')
            roc_curves_120_125[i].GetYaxis().SetTitle('H(120) Rejection')
        else:
            roc_curves_120_125[i].Draw('C')
    roc_legend.Draw('same')
    roc_legend.SetTextSize(0.0375)
    roc_label.Draw('same')
    canvas.SaveAs(_outputFilePath+'/roc_120.png')

#========== H130 vs. H125 ROC ==========#
if _useH:
    roc_curves_125_130=[roccurve.roccurve(H_hists[130][i],H_hists[125][i]) for i in range(3)]
    roc_legend=ROOT.TLegend(0.2,0.2,0.4,0.4)

    roc_curves_125_130[0].SetLineWidth(2)
    roc_curves_125_130[0].SetLineColor(ROOT.kRed)
    roc_legend.AddEntry(roc_curves_125_130[0],'BRT Mass','L')

    roc_curves_125_130[1].SetLineWidth(2)
    roc_curves_125_130[1].SetLineColor(ROOT.kBlack)
    roc_legend.AddEntry(roc_curves_125_130[1],'MMC Mass','L')

    roc_curves_125_130[2].SetLineWidth(2)
    roc_curves_125_130[2].SetLineColor(ROOT.kBlue)
    roc_legend.AddEntry(roc_curves_125_130[2],'Collinear Mass','L')

    for i in range(3):
        if i==0:
            roc_curves_125_130[i].Draw('AC')
            roc_curves_125_130[i].GetXaxis().SetRangeUser(0.,1.)
            roc_curves_125_130[i].GetYaxis().SetRangeUser(0.,1.)
            roc_curves_125_130[i].SetTitle('H(125) vs H(130) ROC Curve')
            roc_curves_125_130[i].GetXaxis().SetTitle('H(130) Acceptance')
            roc_curves_125_130[i].GetYaxis().SetTitle('H(125) Rejection')
        else:
            roc_curves_125_130[i].Draw('C')
    roc_legend.Draw('same')
    roc_legend.SetTextSize(0.0375)
    roc_label.Draw('same')
    canvas.SaveAs(_outputFilePath+'/roc_130.png')

#========== H120 vs. H130 ROC ==========#
if _useH:
    roc_curves_120_130=[roccurve.roccurve(H_hists[130][i],H_hists[120][i]) for i in range(3)]
    roc_legend=ROOT.TLegend(0.2,0.2,0.4,0.4)
    
    roc_curves_120_130[0].SetLineWidth(2)
    roc_curves_120_130[0].SetLineColor(ROOT.kRed)
    roc_legend.AddEntry(roc_curves_120_130[0],'BRT Mass','L')
    
    roc_curves_120_130[1].SetLineWidth(2)
    roc_curves_120_130[1].SetLineColor(ROOT.kBlack)
    roc_legend.AddEntry(roc_curves_120_130[1],'MMC Mass','L')

    roc_curves_120_130[2].SetLineWidth(2)
    roc_curves_120_130[2].SetLineColor(ROOT.kBlue)
    roc_legend.AddEntry(roc_curves_120_130[2],'Collinear Mass','L')
    
    for i in range(3):
        if i==0:
            roc_curves_120_130[i].Draw('AC')
            roc_curves_120_130[i].GetXaxis().SetRangeUser(0.,1.)
            roc_curves_120_130[i].GetYaxis().SetRangeUser(0.,1.)
            roc_curves_120_130[i].SetTitle('H(120) vs H(130) ROC Curve')
            roc_curves_120_130[i].GetXaxis().SetTitle('H(130) Acceptance')
            roc_curves_120_130[i].GetYaxis().SetTitle('H(120) Rejection')
        else:
            roc_curves_120_130[i].Draw('C')
    roc_legend.Draw('same')
    roc_legend.SetTextSize(0.0375)
    roc_label.Draw('same')
    canvas.SaveAs(_outputFilePath+'/roc_130_120.png')

#========== H100 vs. H150 ROC ==========#
if _useH:
    roc_curves_100_150=[roccurve.roccurve(H_hists[150][i],H_hists[100][i]) for i in range(3)]
    roc_legend=ROOT.TLegend(0.2,0.2,0.4,0.4)

    roc_curves_100_150[0].SetLineWidth(2)
    roc_curves_100_150[0].SetLineColor(ROOT.kRed)
    roc_legend.AddEntry(roc_curves_100_150[0],'BRT Mass','L')

    roc_curves_100_150[1].SetLineWidth(2)
    roc_curves_100_150[1].SetLineColor(ROOT.kBlack)
    roc_legend.AddEntry(roc_curves_100_150[1],'MMC Mass','L')

    roc_curves_100_150[2].SetLineWidth(2)
    roc_curves_100_150[2].SetLineColor(ROOT.kBlue)
    roc_legend.AddEntry(roc_curves_100_150[2],'Collinear Mass','L')

    for i in range(3):
        if i==0:
            roc_curves_100_150[i].Draw('AC')
            roc_curves_100_150[i].GetXaxis().SetRangeUser(0.,1.)
            roc_curves_100_150[i].GetYaxis().SetRangeUser(0.,1.)
            roc_curves_100_150[i].SetTitle('H(100) vs H(150) ROC Curve')
            roc_curves_100_150[i].GetXaxis().SetTitle('H(150) Acceptance')
            roc_curves_100_150[i].GetYaxis().SetTitle('H(100) Rejection')
        else:
            roc_curves_100_150[i].Draw('C')
    roc_legend.Draw('same')
    roc_legend.SetTextSize(0.0375)
    roc_label.Draw('same')
    canvas.SaveAs(_outputFilePath+'/roc_150_100.png')
    canvas.Clear()

#========== Variable Plots ==========#
if _variablePlots:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        ks_tests={}
        varlegend=ROOT.TLegend(0.7,0.7,0.9,0.9)
        varlegend.AddEntry(variable_hists_Truth[125]['met_et'],'Truth Level','l')
        varlegend.AddEntry(variable_hists_H[125]['met_et'],'Full Simulation','l')
        for i in _masses:
            for k in variable_hists_H[i].keys():
                stack=ROOT.THStack()
                variable_hists_H[i][k].SetLineColor(ROOT.kRed)
                variable_hists_Truth[i][k].SetLineColor(ROOT.kBlack)
                variable_hists_H[i][k].Sumw2()
                variable_hists_H[i][k].Scale(1./variable_hists_H[i][k].GetEntries())
                variable_hists_Truth[i][k].Sumw2()
                variable_hists_Truth[i][k].Scale(1./variable_hists_Truth[i][k].GetEntries())
                stack.Add(variable_hists_H[i][k])
                stack.Add(variable_hists_Truth[i][k])
                if k in ks_tests.keys():
                    ks_tests[k]+=variable_hists_H[i][k].KolmogorovTest(variable_hists_Truth[i][k])/float(len(_masses))
                else:
                    ks_tests[k]=variable_hists_H[i][k].KolmogorovTest(variable_hists_Truth[i][k])/float(len(_masses))
                stack.Draw('nostack hist')
                varlegend.Draw('same')
                varlegend.SetTextSize(0.0375)
                if i==125 and k=="met_et":
                    met_label=ROOT.TLatex(0.55,0.5,"#splitline{H(125)#rightarrow#tau_{had}#tau_{had} Boosted}{#it{#bf{ATLAS}} Work in Progress}")
                    met_label.SetNDC()
                    met_label.Draw('same')
                    stack.GetXaxis().SetTitle('Missing Transverse Energy (GeV)')
                    stack.GetYaxis().SetTitle('Events / 4 GeV')
                canvas.SaveAs(_outputFilePath+'/variableplots/'+k+'_'+str(i)+'.png')
                canvas.Clear()
        print "P-values for variable comparisons at truth and reco level: "
        f.write("P-values for variable comparisons at truth and reco level: \n")
        for k in sorted(ks_tests, key=ks_tests.get, reverse=True):
            print k,ks_tests[k]
            f.write(k+'    '+str(ks_tests[k]))
        print "\n"
        f.write('\n')

#========== Calculate Metrics ==========#
#P-value
if _useH and _useZ:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        print "Performance Metrics:\n"
        f.write("Performance Metrics:\n")
        chi2=[Z_hists[i].Chi2Test(H_hists[125][i],"UU NORM CHI2") for i in range(3)]
        f.write("Z vs. H(125) Chi^2: ")
        print "Z vs. H(125) Chi^2: "
        f.write("BRT: "+str(chi2[0]))
        print "BRT: ", chi2[0]
        f.write("MMC: "+str(chi2[1]))
        print "MMC: ", chi2[1]
        f.write("Collinear: "+str(chi2[2]))
        print "Collinear: ", chi2[2]
        if max(chi2)==chi2[0]:
            print "Best distinguishing power given by BRT."
            print ""
            f.write("Best distinguishing power given by BRT.\n")
        if max(chi2)==chi2[1]:
            print "Best distinguishing power given by MMC."
            print ""
            f.write("Best distinguishing power given by MMC.\n")
        if max(chi2)==chi2[2]:
            print "Best distinguishing power given by collinear."
            print ""
            f.write("Best distinguishing power given by collinear.\n")

#P-value
if _useH and _useZ:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        hists_temp=[H_hists[125][i].Clone('copy_'+str(i)) for i in range(3)]
        for i in range(3):
            hists_temp[i].Add(Z_hists[i])
            hists_temp[i].SetLineColor(ROOT.kGreen)
            hists_temp[i].Draw('hist')
            Z_hists[i].Draw('same hist')
            if i==0:
                canvas.SaveAs(_outputFilePath+'/Z_H_125_BRT.png')
            elif i==1:
                canvas.SaveAs(_outputFilePath+'/Z_H_125_MMC.png')
            else:
                canvas.SaveAs(_outputFilePath+'/Z_H_125_Col.png')
            canvas.Clear()
        chi2=[Z_hists[i].Chi2Test(hists_temp[i],"WW CHI2/NDF") for i in range(3)]
        f.write("Z + H(125) vs. Z P-value: ")
        print "Z + H(125) vs. Z P-value: "
        f.write("BRT: "+str(chi2[0]))
        print "BRT: ", chi2[0]
        f.write("MMC: "+str(chi2[1]))
        print "MMC: ", chi2[1]
        f.write("Collinear: "+str(chi2[2]))
        print "Collinear: ", chi2[2]
        if min(chi2)==chi2[0]:
            print "Best distinguishing power given by BRT."
            print ""
            f.write("Best distinguishing power given by BRT.\n")
        if min(chi2)==chi2[1]:
            print "Best distinguishing power given by MMC."
            print ""
            f.write("Best distinguishing power given by MMC.\n")
        if min(chi2)==chi2[2]:
            print "Best distinguishing power given by collinear."
            print ""
            f.write("Best distinguishing power given by collinear.\n")    
#Resolution at 125 GeV:
if _useH:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        res_125=[i/125. for i in H_rms[125]]
        print "Resolution @ 125 GeV: "
        f.write("Resolution @ 125 GeV: ")
        print "BRT: ", res_125[0]
        f.write("BRT: "+ str(res_125[0]))
        print "MMC: ", res_125[1]
        f.write("MMC: "+str(res_125[1]))
        print "Collinear: ", res_125[2]
        f.write("Collinear: "+str(res_125[2]))
        if min(res_125)==res_125[0]:
            print "Best distinguishing power given by BRT."
            print ""
            f.write("Best distinguishing power given by BRT.")
        if min(res_125)==res_125[1]:
            print "Best distinguishing power given by MMC."
            print ""
            f.write("Best distinguishing power given by MMC.\n")
        if min(res_125)==res_125[2]:
            print "Best distinguishing power given by collinear."
            print ""
            f.write("Best distinguishing power given by collinear.\n")
print Z_mean
print H_mean[125]
#Resolution at Z mass:
if _useZ:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        res_Z=[i/_Zmass for i in Z_rms]
        print "Resolution @ Z mass: "
        f.write("Resolution @ Z mass: ")
        print "BRT: ", res_Z[0]
        f.write("BRT: "+str(res_Z[0]))
        print "MMC: ", res_Z[1]
        f.write("MMC: "+str(res_Z[1]))
        print "Collinear: ", res_Z[2]
        f.write("Collinear: "+str(res_Z[2]))
        if min(res_Z)==res_Z[0]:
            print "Best distinguishing power given by BRT."
            print ""
            f.write("Best distinguishing power given by BRT.\n")
        if min(res_Z)==res_Z[1]:
            print "Best distinguishing power given by MMC."
            print ""
            f.write("Best distinguishing power given by MMC.\n")
        if min(res_Z)==res_Z[2]:
            print "Best distinguishing power given by collinear."
            print ""
            f.write("Best distinguishing power given by collinear.\n")

#Z vs. H(125) slope
if _useZ and _useH:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        slope=[(H_mean[125][i]-Z_mean[i])/(125.-_Zmass) for i in range(3)]
        print "Slope between Z and H(125): "
        f.write("Slope between Z and H(125): ")
        print "BRT: ", slope[0]
        f.write("BRT: "+str(slope[0]))
        print "MMC: ", slope[1]
        f.write("MMC: "+str(slope[1]))
        print "Collinear: ", slope[2]
        f.write("Collinear: "+str(slope[2]))
        if max(slope)==slope[0]:
            print "Best distinguishing power given by BRT."
            print ""
            f.write("Best distinguishing power given by BRT.\n")
        if max(slope)==slope[1]:
            print "Best distinguishing power given by MMC."
            print ""
            f.write("Best distinguishing power given by MMC.\n")
        if max(slope)==slope[2]:
            print "Best distinguishing power given by collinear."
            print ""
            f.write("Best distinguishing power given by collinear.\n")

#Fitted slope
if _useZ and _useH:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        slope_fitted=[]
        for i in range(3):
            meanGraphs[i].Fit('pol1')
            fitfunc=meanGraphs[i].GetFunction('pol1')
            slope_fitted.append(fitfunc.GetParameter(1))
        print "Slope between Z and H(125): "
        f.write("Slope between Z and H(125): ")
        print "BRT: ", slope_fitted[0]
        f.write("BRT: "+str(slope_fitted[0]))
        print "MMC: ", slope_fitted[1]
        f.write("MMC: "+str(slope_fitted[1]))
        print "Collinear: ", slope_fitted[2]
        f.write("Collinear: "+str(slope_fitted[2]))
        if max(slope_fitted)==slope_fitted[0]:
            print "Best distinguishing power given by BRT."
            print ""
            f.write("Best distinguishing power given by BRT.\n")
        if max(slope_fitted)==slope_fitted[1]:
            print "Best distinguishing power given by MMC."
            print ""
            f.write("Best distinguishing power given by MMC.\n")
        if max(slope_fitted)==slope_fitted[2]:
            print "Best distinguishing power given by collinear."
            print ""
            f.write("Best distinguishing power given by collinear.\n")

#Slope/Resolution at 125 GeV:
if _useH and _useZ:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        slope_res_125=[((H_mean[125][i]-Z_mean[i])/(125.-_Zmass))/(H_rms[125][i]/125.) for i in range(3)]
        print "Slope/Resolution @ 125 GeV: "
        f.write("Slope/Resolution @ 125 GeV: ")
        print "BRT: ", slope_res_125[0]
        f.write("BRT: "+str(slope_res_125[0]))
        print "MMC: ", slope_res_125[1]
        f.write("MMC: "+str(slope_res_125[1]))
        print "Collinear: ", slope_res_125[2]
        f.write("Collinear: "+str(slope_res_125[2]))
        if max(slope_res_125)==slope_res_125[0]:
            print "Best distinguishing power given by BRT."
            print ""
            f.write("Best distinguishing power given by BRT.\n")
        if max(slope_res_125)==slope_res_125[1]:
            print "Best distinguishing power given by MMC."
            print ""
            f.write("Best distinguishing power given by MMC.\n")
        if max(slope_res_125)==slope_res_125[2]:
            print "Best distinguishing power given by collinear."
            print ""
            f.write("Best distinguishing power given by collinear.\n")

#Slope/Resolution at Zmass:
if _useH and _useZ:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        slope_res_Z=[((H_mean[125][i]-Z_mean[i])/(125.-_Zmass))/(Z_rms[i]/_Zmass) for i in range(3)]
        print "Slope/Resolution @ Z mass: "
        f.write("Slope/Resolution @ Z mass: ")
        print "BRT: ", slope_res_Z[0]
        f.write("BRT: "+str(slope_res_Z[0]))
        print "MMC: ", slope_res_Z[1]
        f.write("MMC: "+str(slope_res_Z[1]))
        print "Collinear: ", slope_res_Z[2]
        f.write("Collinear: "+str(slope_res_Z[2]))
        if max(slope_res_Z)==slope_res_Z[0]:
            print "Best distinguishing power given by BRT."
            print ""
            f.write("Best distinguishing power given by BRT.\n")
        if max(slope_res_Z)==slope_res_Z[1]:
            print "Best distinguishing power given by MMC."
            print ""
            f.write("Best distinguishing power given by MMC.\n")
        if max(slope_res_Z)==slope_res_Z[2]:
            print "Best distinguishing power given by collinear."
            print ""
            f.write("Best distinguishing power given by collinear.\n")

#ROC Integral
if _useH and _useZ:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        integrals_Z_125=[i.Integral() for i in roc_curves_Z_125]
        print "Z vs. H(125) ROC curve areas: "
        f.write("Z vs. H(125) ROC curve areas: ")
        print "BRT: ", integrals_Z_125[0]
        f.write("BRT: "+str(integrals_Z_125[0]))
        print "MMC: ", integrals_Z_125[1]
        f.write("MMC: "+str(integrals_Z_125[1]))
        print "Collinear: ", integrals_Z_125[2]
        f.write("Collinear: "+str(integrals_Z_125[2]))
        if max(integrals_Z_125)==integrals_Z_125[0]:
            print "Best distinguishing power given by BRT."
            print ""
            f.write("Best distinguishing power given by BRT.\n")
        if max(integrals_Z_125)==integrals_Z_125[1]:
            print "Best distinguishing power given by MMC."
            print ""
            f.write("Best distinguishing power given by MMC.\n")
        if max(integrals_Z_125)==integrals_Z_125[2]:
            print "Best distinguishing power given by collinear."
            print ""
            f.write("Best distinguishing power given by collinear.\n")

#Signal Acceptance @ Signal Acceptance=Background Rejection
if _useH and _useZ:
    with open(_outputFilePath+'/output.txt', 'a') as f:
        cutoff_error=0.0001
        SA=[]
        for i in range(3):
            error=1.
            x=0.8
            counter=0
            while error>cutoff_error:
                y=roc_curves_Z_125[i].Eval(x)
                error=abs(y-x)
                x=y
                counter=counter+1
                if counter>1000:
                    break
            if counter>1000:
                f.write("Insufficient statistics for SA=BR intercept.\n")
                print "Insufficient statistics for SA=BR intercept.\n"
                break
            SA.append(x)
        if counter<=1000:
            print "Z vs. H(125) ROC curve SA=BR intercept: "
            f.write("Z vs. H(125) ROC curve SA=BR intercept: ")
            print "BRT: ", SA[0]
            f.write("BRT: "+str(SA[0]))
            print "MMC: ", SA[1]
            f.write("MMC: "+str(SA[1]))
            print "Collinear: ", SA[2]
            f.write("Collinear: "+str(SA[2]))
            if max(SA)==SA[0]:
                print "Best distinguishing power given by BRT."
                print ""
                f.write("Best distinguishing power given by BRT.\n")
            if max(SA)==SA[1]:
                print "Best distinguishing power given by MMC."
                print ""
                f.write("Best distinguishing power given by MMC.\n")
            if max(SA)==SA[2]:
                print "Best distinguishing power given by collinear."
                print ""
                f.write("Best distinguishing power given by collinear.\n")

#========== Write hists to file ==========#
outputfile=ROOT.TFile(_outputFilePath+'/graph.root','RECREATE')
for i in range(3):
    if _useZ:
        Z_hists[i].Write()
    if _useH:
        for j in _masses:
            H_hists[j][i].Write()
        if _useZ:
            roc_curves_Z_125[i].Write()
        roc_curves_120_125[i].Write()
        roc_curves_125_130[i].Write()
        roc_curves_100_150[i].Write()
if _variablePlots:
    for imass in _masses:
        for k in variable_hists_Truth[imass].keys():
            variable_hists_Truth[imass][k].Write()
            variable_hists_H[imass][k].Write()
outputfile.Close()
