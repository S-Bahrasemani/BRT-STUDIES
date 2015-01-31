####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
## Code: TrainingAndTesting.py                                                                     ##
## Authors: Natascha Hedrich (hedrich.natascha(AT)gmail.com)                                       ##
##          Andres Tanasijczuk (andres.tanasijczuk(AT)gmail.com)                                   ##
##          Matthew Quenneville (mattq1(AT)gmail.com)                                              ##
## Purpose: Use of boosted regression trees to reconstruct the invariant mass in X->tautau decays. ##
##          The training is done on MC@NLO Higgs ggH samples, with masses between 40 and 200 GeV.  ##
##          The testing is done also on these Higgs samples, but on orthogonal events of course.   ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################

from ROOT import *
import os, sys, math
import pickle
import array
from random import gauss
from random import uniform
from datetime import datetime
from time import sleep
from threading import Timer
import thread
import select, subprocess
import collinearmass

TMVA.Tools.Instance()

#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
## Auxiliary functions for logical AND, OR and XOR (exclusive OR).                                 ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################

def AND(*args):
    for i in range(len(args)):
        if not args[i]: return False
    return True

def OR(*args):
    for i in range(len(args)):
        if args[i]: return True
    return False

def XOR(*args):
    passed = False
    for i in range(len(args)):
        if args[i] and passed: return False
        if args[i] and not passed: passed = True
    return passed

#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
##                                          SETTINGS                                               ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################

#The _useHtt Boolean has been removed, as I do not see the purpose in training without Htt events, and have thus neglected to keep this option up to date.
_useZtt = False
#I am unsure of where this value was obtained, but it is not the correct Z mass, and I have thus replaced it.
#_Zmass = 91.55
_Zmass=91.1876

_Zmass_window = [_Zmass-3.738,_Zmass+3.738]

#Make plots using average mass instead of nominal mass. This should generally be set to True, as this seems to correct biases caused by varying distribution shapes.
_useAverageMass=True

#ditau_m, alpha, alpha_vis
_target='ditau_m'

_fillTrainTestTrees=False
_doTraining =True
_doTesting  = True

if _doTesting:
    _testWOCalibration = True
    _doCalibration     = False
    _testCalibration   = False
    if _doCalibration and not(_target=='ditau_m'):
        print 'Calibration is only available for regression on ditau mass.'
        sys.exit()

    if not OR(_testWOCalibration,_doCalibration,_testCalibration):
        print 'Error: You set "_doTesting" = True, so you also need to set "_testWOCalibration", "_doCalibration" or "_testCalibration" = True.'
        sys.exit()

_useBosonVariables           = False
_useTauVariables             = True
_useTauDecayVariables        = False
_useAnalysisVariables        = False
_useSmearedAnalysisVariables = False
#To add
_useFullSimVariables         = False

if not XOR(_useBosonVariables,_useTauVariables,_useTauDecayVariables,_useAnalysisVariables,_useSmearedAnalysisVariables,_useFullSimVariables):
    print 'Error: Only one of the variables type can/must be True.'
    sys.exit()
if not  OR(_useBosonVariables,_useTauVariables,_useTauDecayVariables,_useAnalysisVariables,_useSmearedAnalysisVariables,_useFullSimVariables):
    print 'Error: You need to set one of the variables type.'
    sys.exit()

_doLepLep = False
_doLepHad = False
_doHadHad = True

if OR(_doTraining,_doTesting):
    if OR(_doLepLep,_doLepHad,_doHadHad) and not XOR(_doLepLep,_doLepHad,_doHadHad):
        print 'Error: Only one of the tau decay channels can be True.'
        sys.exit()
    if OR(_useTauDecayVariables,_useAnalysisVariables,_useSmearedAnalysisVariables) and not OR(_doLepLep,_doLepHad,_doHadHad):
        print 'Error: You need to set the di-tau decay channel.'
        sys.exit()

_applyPreselCuts = False
_useAnalysisPresel=False
if OR(_doTraining,_doTesting):
    if _applyPreselCuts:
        if not OR(_useAnalysisVariables,_useSmearedAnalysisVariables):
            print 'Error: preselection cuts can be applied only to analysis type variables.'
            sys.exit()

if _useAnalysisPresel and not(_applyPreselCuts):
    print 'Error: Analysis preselection can only be applied alongside the applyPreselCuts option'
    sys.exit()
## Cuts
##--------------------
## |eta_lepton| < 2.4
## pT_lepton > 15 GeV
## MET > 20 GeV
##------------------
##Analysis Preselection
##|eta_lepton|<2.5
##pT_had_1>35000
##pT_had_2>25000
##MET >20000

_mass_points_train = range(41,201)

_mass_points_calibration = [41,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200]

_mass_points_test = range(41,201)

if _doTesting and _testCalibration:
    for m in _mass_points_calibration:
        if m in _mass_points_test:
            _mass_points_test.pop(_mass_points_test.index(m))
_mass_points_test.sort()

_mass_points_all = []
for m in _mass_points_train + _mass_points_calibration + _mass_points_test:
    if m not in _mass_points_all:
        _mass_points_all.append(m)
_mass_points_all.sort()

_HttSamplesDirName = '/cluster/data04/mquennev/higgs/rootfiles'
_ZttSamplesDirName = _HttSamplesDirName
_TreeNameInput     = 'Tree'

_treesDirName   = '/cluster/data04/mquennev/higgs/trees'
_testDirName    = '/cluster/data04/mquennev/higgs/test'
_weightsDirName = '/cluster/data04/mquennev/higgs/weights'
_TreeNameTrain  = 'TreeTrain'
_TreeNameTest   = 'TreeTest'

variables = {}
if _useBosonVariables:
    variables['boson_pt']  = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['boson_eta'] = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['boson_phi'] = ['rad','F', -3.15, 3.15,  array.array('f',[0]),TBranch(),  -3.15, 3.15,   60]
    variables['boson_E']   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
if _useTauVariables:
    variables['taup_pt']   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['taup_eta']  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['taup_phi']  = ['rad','F', -3.15, 3.15,  array.array('f',[0]),TBranch(),  -3.15, 3.15,   60]
    variables['taup_E']    = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['taum_pt']   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['taum_eta']  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['taum_phi']  = ['rad','F', -3.15, 3.15,  array.array('f',[0]),TBranch(),  -3.15, 3.15,   60]
    variables['taum_E']    = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
if _useTauDecayVariables:
    variables['tau1_decay_nutau_pt']   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['tau1_decay_nutau_eta']  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['tau1_decay_nutau_phi']  = ['rad','F', -3.15, 3.15,  array.array('f',[0]),TBranch(),  -3.15, 3.15,   60]
    variables['tau1_decay_nutau_E']    = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['tau2_decay_nutau_pt']   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['tau2_decay_nutau_eta']  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['tau2_decay_nutau_phi']  = ['rad','F', -3.15, 3.15,  array.array('f',[0]),TBranch(),  -3.15, 3.15,   60]
    variables['tau2_decay_nutau_E']    = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['tau1_decay_lep_pt']     = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['tau1_decay_lep_eta']    = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['tau1_decay_lep_phi']    = ['rad','F', -3.15, 3.15,  array.array('f',[0]),TBranch(),  -3.15, 3.15,   60]
    variables['tau1_decay_lep_E']      = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['tau2_decay_lep_pt']     = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['tau2_decay_lep_eta']    = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['tau2_decay_lep_phi']    = ['rad','F', -3.15, 3.15,  array.array('f',[0]),TBranch(),  -3.15, 3.15,   60]
    variables['tau2_decay_lep_E']      = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    if _doLepLep:
        variables['tau1_decay_nulep_pt']   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
        variables['tau1_decay_nulep_eta']  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
        variables['tau1_decay_nulep_phi']  = ['rad','F', -3.15, 3.15,  array.array('f',[0]),TBranch(),  -3.15, 3.15,   60]
        variables['tau1_decay_nulep_E']    = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
        variables['tau2_decay_nulep_pt']   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
        variables['tau2_decay_nulep_eta']  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
        variables['tau2_decay_nulep_phi']  = ['rad','F', -3.15, 3.15,  array.array('f',[0]),TBranch(),  -3.15, 3.15,   60]
        variables['tau2_decay_nulep_E']    = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    if _doLepHad:
        variables['tau1_decay_nulep_pt']   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
        variables['tau1_decay_nulep_eta']  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
        variables['tau1_decay_nulep_phi']  = ['rad','F', -3.15, 3.15,  array.array('f',[0]),TBranch(),  -3.15, 3.15,   60]
        variables['tau1_decay_nulep_E']    = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
if _useAnalysisVariables:
    variables['lep1_pt']                      = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['lep1_eta']                     = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['lep2_pt']                      = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['lep2_eta']                     = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['met_et']                       = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['transverse_mass_lep1_lep2']    = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 50]
    variables['transverse_mass_lep1_met']     = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 50]
    variables['transverse_mass_lep2_met']     = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 50]
    variables['dphi_lep1_met']                = ['rad','F',  0.00, 3.15,  array.array('f',[0]),TBranch(),   0.00, 3.15,   30]
    variables['dphi_lep2_met']                = ['rad','F',  0.00, 3.15,  array.array('f',[0]),TBranch(),   0.00, 3.15,   30]
    variables['dphi_lep_lep']                 = ['rad','F',  0.00, 3.15,  array.array('f',[0]),TBranch(),   0.00, 3.15,   30]
    variables['deta_lep_lep']                 = ['',   'F',  0.00,20.00,  array.array('f',[0]),TBranch(),   0.00,20.00,  200]
    variables['dR_lep_lep']                   = ['',   'F',  0.00,25.00,  array.array('f',[0]),TBranch(),   0.00,20.00,  200]
    variables['ptsum_lep1_lep2_met']          = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,300000., 30]
    variables['ptsum_lep1_lep2']              = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 25]
    variables['pttot_lep1_lep2_met']          = ['',   'F',  0.00, 2.00,  array.array('f',[0]),TBranch(),   0.00, 1.10,   22]
    variables['pttot_lep1_lep2']              = ['',   'F',  0.00, 2.00,  array.array('f',[0]),TBranch(),   0.00, 1.10,   22]
    variables['ptdiff_lep1_lep2']             = ['',   'F',  0.00, 2.00,  array.array('f',[0]),TBranch(),   0.00, 1.10,   22]

if _useSmearedAnalysisVariables:
    variables['lep1_pt_sm']                   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['lep1_eta_sm']                  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['lep2_pt_sm']                   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['lep2_eta_sm']                  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
    variables['met_et_sm']                    = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
    variables['transverse_mass_lep1_lep2_sm'] = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 50]
    variables['transverse_mass_lep1_met_sm']  = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 50]
    variables['transverse_mass_lep2_met_sm']  = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 50]
    variables['dphi_lep1_met_sm']             = ['rad','F',  0.00, 3.15,  array.array('f',[0]),TBranch(),   0.00, 3.15,   30]
    variables['dphi_lep2_met_sm']             = ['rad','F',  0.00, 3.15,  array.array('f',[0]),TBranch(),   0.00, 3.15,   30]
    variables['dphi_lep_lep_sm']              = ['rad','F',  0.00, 3.15,  array.array('f',[0]),TBranch(),   0.00, 3.15,   30]
    variables['deta_lep_lep_sm']              = ['',   'F',  0.00,20.00,  array.array('f',[0]),TBranch(),   0.00,20.00,  200]
    variables['dR_lep_lep_sm']                = ['',   'F',  0.00,25.00,  array.array('f',[0]),TBranch(),   0.00,20.00,  200]
    variables['ptsum_lep1_lep2_met_sm']       = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,300000., 30]
    variables['ptsum_lep1_lep2_sm']           = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 25]
    variables['pttot_lep1_lep2_met_sm']       = ['',   'F',  0.00, 2.00,  array.array('f',[0]),TBranch(),   0.00, 1.10,   22]
    variables['pttot_lep1_lep2_sm']           = ['',   'F',  0.00, 2.00,  array.array('f',[0]),TBranch(),   0.00, 1.10,   22]
    variables['ptdiff_lep1_lep2_sm']          = ['',   'F',  0.00, 2.00,  array.array('f',[0]),TBranch(),   0.00, 1.10,   22]

channel_string = ''
if _doLepLep: channel_string = 'H #rightarrow #tau_{lep}#tau_{lep}'
if _doLepHad: channel_string = 'H #rightarrow #tau_{lep}#tau_{had}'
if _doHadHad: channel_string = 'H #rightarrow #tau_{had}#tau_{had}'


#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
## Function to make the files with training and testing trees.                                     ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################


def FillTrainTestTrees(iFileName,oFileName,nTrainEvtsMax):

    print "<--  input file: "+iFileName
    print "--> output file: "+oFileName
    iFile = TFile.Open(iFileName)
    iTree = iFile.Get(_TreeNameInput)
    oDirName = oFileName[::-1].split("/",1)[-1][::-1]
    os.system("mkdir -p "+oDirName)
    oFile = TFile(oFileName,"RECREATE")
    oTreeTrain = iTree.CloneTree(0)
    oTreeTrain.SetName(_TreeNameTrain)
    oTreeTest  = iTree.CloneTree(0)
    oTreeTest.SetName(_TreeNameTest)
    nTrainEvts = 0
    nTestEvts  = 0
    nEntries = iTree.GetEntries()
    nTrainEvtsMax = round(nEntries/2)
    print "number of entries in original tree = %s" %nEntries

    for ientry in xrange(nEntries):
        ## Get the next tree in the chain and verify.
        if iTree.LoadTree(ientry) <  0: break
        ## Copy next entry into memory and verify.
        if iTree.GetEntry(ientry) <= 0: continue
        if _useBosonVariables or _useTauVariables:
            if nTrainEvts < nTrainEvtsMax/10.:
                oTreeTrain.Fill()
                nTrainEvts += 1
            else:
                oTreeTest.Fill()
                nTestEvts += 1
        else:
            if nTrainEvts < nTrainEvtsMax:
                oTreeTrain.Fill()
                nTrainEvts += 1
            else:
                oTreeTest.Fill()
                nTestEvts  += 1
    oTreeTrain.AutoSave()
    oTreeTest.AutoSave()
    iFile.Close()
    oFile.Close()
    print "Saved %i events in the train tree" %nTrainEvts
    print "Saved %i events in the test tree" %nTestEvts


#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
##                                        BRT TRAINING                                             ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################


def Training(factoryName,methodName,trainParams):

    ## Group all the training samples.
    TrainSamples = []
    if _useZtt:
        TrainSamples.append(ZttSample())
    for mass in _mass_points_train:
        if _useZtt and _Zmass_window[0] < mass < _Zmass_window[-1]: continue
        TrainSamples.append(HttSample(mass))

    ## Add the training samples to a chain to be used by the factory.
    TrainTree = TChain(_TreeNameTrain)
    for sample in TrainSamples:
        TrainTree.Add(sample)
    entries=TrainTree.GetEntries()

    ## Create the output file where all variables will be stored.
    training_parameters = '!H:V:BoostType=AdaBoostR2:SeparationType=RegressionVariance'
    for key, value in trainParams.iteritems():
        training_parameters += ':'+str(key)+'='+str(value)
        methodName += '_'+str(key)+str(value)
    outSubDirName = factoryName+'_'+methodName
    outdir = GetDir(_testDirName)+'/'+outSubDirName
    os.system('mkdir -p '+outdir)
    outFile = TFile(outdir+'/TrainingOutput.root','RECREATE')

    ## Set the directory where to save the weights file.
    (TMVA.gConfig().GetIONames()).fWeightFileDir = GetDir(_weightsDirName)
    os.system('mkdir -p '+GetDir(_weightsDirName))

    ## Initialize the TMVA factory.
    factory = TMVA.Factory(factoryName,outFile,'V:!Silent:Color') #:DrawProgressBar')

    ## Add variables to the factory.
    for varName, var in sorted(variables.iteritems()):
        factory.AddVariable(varName,varName,var[0],var[1],var[2],var[3])
    #factory.AddSpectator("collinear_mass")
    #factory.AddSpectator("visible_mass")

    ## Add the target.
    factory.AddTarget(_target)

    #factory.SetWeightExpression('array3')

    ## Add the input files to the factory.
    factory.AddRegressionTree(TrainTree)    

    factory.PrepareTrainingAndTestTree(TCut(Cuts(True,False)),'nTrain_Regression=0:nTest_Regression=1:SplitMode=Random:NormMode=NumEvents:!V')
    factory.BookMethod(TMVA.Types.kBDT,methodName,training_parameters)

    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()

    outFile.Close()
    gROOT.LoadMacro("$ROOTSYS/tmva/test/TMVAGui.C")
    TMVAGui(outdir+"/TrainingOutput.root")
    raw_input('Press Enter to exit')
    ## Make a plot with the true mass in the training samples.
    canvas = TCanvas('canvas','canvas',500,525)
    canvas.GetPad(0).SetLeftMargin(0.15)
    canvas.GetPad(0).SetRightMargin(0.05)
    hists = {}
    mass_min = 40
    mass_max = 200
    hists['mass_true_all_train'] = TH1F('h_mass_true_all_train','m_{true} for all samples in train tree',(mass_max-mass_min+19)*10,mass_min-9.5,mass_max+9.5)
    hists['mass_true_all_train'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['mass_true_all_train'].GetYaxis().SetTitle('Events')
    hists['mass_true_all_train'].GetXaxis().SetTitleOffset(1.3)
    hists['mass_true_all_train'].GetYaxis().SetTitleOffset(2.0)
    for sample in TrainSamples:
        file = TFile.Open(sample)
        tree = file.Get(_TreeNameTrain)
        nEntries = tree.GetEntriesFast()
        for ientry in xrange(nEntries):
            ## Get the next tree in the chain and verify
            if tree.LoadTree(ientry) <  0: break
            ## Copy next entry into memory and verify
            if tree.GetEntry(ientry) <= 0: continue
            ## Apply decay channel and selection cuts
            if not PassedCuts(False,True,tree): continue
            mass_true = tree.ditau_m/1000.
            hists['mass_true_all_train'].Fill(mass_true)
        file.Close()
    hists['mass_true_all_train'].Draw()
    canvas.SaveAs(outdir+'/'+hists['mass_true_all_train'].GetName()+'.png')
    canvas.Clear()


#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
##                                        BRT TESTING                                              ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################


def Testing(xmlFileName,testWOCalibration,doCalibration,testCalibration):

    if testWOCalibration:
        doCalibration     = False
        testCalibration   = False
    if doCalibration:
        testWOCalibration = False
        testCalibration   = False
    if testCalibration:
        testWOCalibration = False
        doCalibration     = False

    outSubDirName = xmlFileName.replace('.weights.xml','')
    outdir = GetDir(_testDirName)+'/'+outSubDirName

    ## Ask user if sure about running the testing
    if testWOCalibration:
         print 'Will test using weights file '+GetDir(_weightsDirName)+'/'+xmlFileName
    elif doCalibration:
        print 'Will do calibration using weights file '+GetDir(_weightsDirName)+'/'+xmlFileName
    elif testCalibration:
        print 'Will test calibration using weights file '+GetDir(_weightsDirName)+'/'+xmlFileName

    os.system('mkdir -p '+outdir)

    if doCalibration:
        calibFile = open(outdir+'/calibration.txt','w')
        calibration = []
    elif testCalibration:
        calibFile = open(outdir+'/calibration.txt','r')
        calibration = pickle.load(calibFile)

    if doCalibration:
        mass_points = _mass_points_calibration
    elif testWOCalibration or testCalibration:
        mass_points = _mass_points_test
    else:
        mass_points = _mass_points_all

    canvas = TCanvas('canvas','canvas',500,525)
    canvas.GetPad(0).SetLeftMargin(0.15)
    canvas.GetPad(0).SetRightMargin(0.05)
    latex = TLatex()
    latex.SetTextFont(42)
    latex.SetTextSize(0.035)
    line = TLine()
    box = TBox()

    files = {}
    trees = {}
    hists = {}
    for mass in mass_points:
        files[str(mass)] = TFile.Open(HttSample(mass))
        trees[str(mass)] = files[str(mass)].Get(_TreeNameTest)
    if _useZtt:
        files["Z"] = TFile.Open(ZttSample())
        trees["Z"] = files["Z"].Get(_TreeNameTest)

    outFile = TFile(outdir+'/TestingOutput.root','RECREATE')

    if testWOCalibration or doCalibration:
        auxNameStr = '_uncalib'
        auxTitleStr = 'uncalibrated '
        mBRTleg = 'm^{uncalib}_{BRT}'
    elif testCalibration:
        auxNameStr = '_calib'
        auxTitleStr = 'calibrated '
        mBRTleg = 'm^{calib}_{BRT}'
    else:
        auxNameStr = ''
        auxTitleStr = ''
        mBRTleg = 'm_{BRT}'

    ## Define the histograms.
    mass_min = 40
    mass_max = 200
    #-----------
    for mass in mass_points:
        hists['mass_BRT_'+str(mass)] = TH1F('h_mass_BRT'+auxNameStr+'_'+str(mass),auxTitleStr+'BRT mass for '+str(mass)+' GeV Higgs sample',(mass_max-mass_min+20),mass_min-10,mass_max+10)
        hists['mass_BRT_'+str(mass)].GetXaxis().SetTitle(mBRTleg+' [GeV]')
        hists['mass_BRT_'+str(mass)].GetYaxis().SetTitle('Events')
        hists['mass_BRT_'+str(mass)].GetXaxis().SetTitleOffset(1.3)
        hists['mass_BRT_'+str(mass)].GetYaxis().SetTitleOffset(2.0)
    #-----------
    if _useZtt:
        hists['mass_BRT_Z'] = TH1F('h_mass_BRT'+auxNameStr+'_Z',auxTitleStr+'BRT mass for Z sample',(mass_max-mass_min+20),mass_min,mass_max)
        hists['mass_BRT_Z'].GetXaxis().SetTitle(mBRTleg+' [GeV]')
        hists['mass_BRT_Z'].GetYaxis().SetTitle('Events')
        hists['mass_BRT_Z'].GetXaxis().SetTitleOffset(1.3)
        hists['mass_BRT_Z'].GetYaxis().SetTitleOffset(2.0)
        hists['mass_true_Z'] = TH1F('h_mass_true_Z'+auxNameStr,'m_{true} for Z sample in test tree',80,70,110)
        hists['mass_true_Z'].GetXaxis().SetTitle('m_{true} [GeV]')
        hists['mass_true_Z'].GetYaxis().SetTitle('Events')
        hists['mass_true_Z'].GetXaxis().SetTitleOffset(1.3)
        hists['mass_true_Z'].GetYaxis().SetTitleOffset(2.0)
    #-----------
    hists['mass_BRT_all'] = TH1F('h_mass_BRT'+auxNameStr+'_all','',(mass_max-mass_min+20),mass_min-10,mass_max+10)
    hists['mass_BRT_all'].GetXaxis().SetTitle(mBRTleg+' [GeV]')
    hists['mass_BRT_all'].GetYaxis().SetTitle('Events')
    hists['mass_BRT_all'].GetXaxis().SetTitleOffset(1.3)
    hists['mass_BRT_all'].GetYaxis().SetTitleOffset(2.0)
    hists['mass_true_all'] = TH1F('h_mass_true_all'+auxNameStr,'m_{true} for all Higgs samples in test tree',(mass_max-mass_min+19)*10,mass_min-9.5,mass_max+9.5)
    hists['mass_true_all'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['mass_true_all'].GetYaxis().SetTitle('Events')
    hists['mass_true_all'].GetXaxis().SetTitleOffset(1.3)
    hists['mass_true_all'].GetYaxis().SetTitleOffset(2.0)
    #-----------
    hists['mass_BRT_bias'] = TH1F('h_bias_mass_BRT'+auxNameStr,'',200,-100,100)
    hists['mass_BRT_bias'].GetXaxis().SetTitle(mBRTleg+' - m_{true} [GeV]')
    hists['mass_BRT_bias'].GetYaxis().SetTitle('Events')
    hists['mass_BRT_bias'].GetXaxis().SetTitleOffset(1.3)
    hists['mass_BRT_bias'].GetYaxis().SetTitleOffset(2.0)
    #-----------
    hists['mean_mass_BRT_vs_mass_true'] = TGraphErrors()
    hists['mean_mass_BRT_vs_mass_true'].SetName('g_mean_mass_BRT'+auxNameStr+'_vs_mass_true')
    hists['mean_mass_BRT_vs_mass_true'].SetMarkerSize(0.3)
    hists['mean_mass_BRT_vs_mass_true'].SetMarkerStyle(20)
    if doCalibration and _polyfit:
        hists['mass_true_vs_mean_mass_BRT'] = TGraphErrors()
        hists['mass_true_vs_mean_mass_BRT'].SetName('g_mass_true'+auxNameStr+'_vs_mean_mass_BRT')
        hists['mass_true_vs_mean_mass_BRT'].SetMarkerSize(0.3)
        hists['mass_true_vs_mean_mass_BRT'].SetMarkerStyle(20)
    #-----------
    hists['max_mass_BRT_vs_mass_true'] = TGraphErrors()
    hists['max_mass_BRT_vs_mass_true'].SetName('g_max_mass_BRT'+auxNameStr+'_vs_mass_true')
    hists['max_mass_BRT_vs_mass_true'].SetMarkerSize(0.3)
    hists['max_mass_BRT_vs_mass_true'].SetMarkerStyle(20)
    #-----------
    hists['med_mass_BRT_vs_mass_true'] = TGraphErrors()
    hists['med_mass_BRT_vs_mass_true'].SetName('g_med_mass_BRT'+auxNameStr+'_vs_mass_true')
    hists['med_mass_BRT_vs_mass_true'].SetMarkerSize(0.3)
    hists['med_mass_BRT_vs_mass_true'].SetMarkerStyle(20)
    #-----------
    hists['rms_mass_BRT_vs_mass_true'] = TGraphErrors()
    hists['rms_mass_BRT_vs_mass_true'].SetName('g_rms_mass_BRT'+auxNameStr+'_vs_mass_true')
    hists['rms_mass_BRT_vs_mass_true'].SetMarkerSize(0.3)
    hists['rms_mass_BRT_vs_mass_true'].SetMarkerStyle(20)
    #-----------
    hists['res_mass_BRT_vs_mass_true'] = TGraphErrors()
    hists['res_mass_BRT_vs_mass_true'].SetName('g_res_mass_BRT'+auxNameStr+'_vs_mass_true')
    hists['res_mass_BRT_vs_mass_true'].SetMarkerSize(0.3)
    hists['res_mass_BRT_vs_mass_true'].SetMarkerStyle(20)
    #-----------
    if testCalibration:
        hists['calib_efficiency_vs_mass_true'] = TGraphErrors()
        hists['calib_efficiency_vs_mass_true'].SetName('g_calib_efficiency_vs_mass_true')
        hists['calib_efficiency_vs_mass_true'].SetMarkerSize(0.3)
        hists['calib_efficiency_vs_mass_true'].SetMarkerStyle(20)
    #-----------
    for varName, var in sorted(variables.iteritems()):
        hists[varName] = TH1F('h_'+varName+auxNameStr,'',var[8],var[6],var[7])
        hists[varName].GetXaxis().SetTitle(varName)
        hists[varName].GetYaxis().SetTitle('Events')
        hists[varName].GetXaxis().SetTitleOffset(1.3)
        hists[varName].GetYaxis().SetTitleOffset(2.0)
    #-----------
    hists['sigma_left_vs_mass_true'] = TGraphErrors()
    hists['sigma_left_vs_mass_true'].SetName('sigma_left'+auxNameStr+'_vs_mass_true')
    hists['sigma_left_vs_mass_true'].SetMarkerSize(0.3)
    hists['sigma_left_vs_mass_true'].SetMarkerStyle(20)
    hists['sigma_left_vs_mass_true'].SetMarkerColor(kRed)
    hists['sigma_left_vs_mass_true'].SetLineColor(kRed)
    hists['sigma_right_vs_mass_true'] = TGraphErrors()
    hists['sigma_right_vs_mass_true'].SetName('sigma_right'+auxNameStr+'_vs_mass_true')
    hists['sigma_right_vs_mass_true'].SetMarkerSize(0.3)
    hists['sigma_right_vs_mass_true'].SetMarkerStyle(20)
    hists['sigma_right_vs_mass_true'].SetMarkerColor(kBlue)
    hists['sigma_right_vs_mass_true'].SetLineColor(kBlue)
    #-----------
    hists['res_left_vs_mass_true'] = TGraphErrors()
    hists['res_left_vs_mass_true'].SetName('res_left'+auxNameStr+'_vs_mass_true')
    hists['res_left_vs_mass_true'].SetMarkerSize(0.3)
    hists['res_left_vs_mass_true'].SetMarkerStyle(20)
    hists['res_left_vs_mass_true'].SetMarkerColor(kRed)
    hists['res_left_vs_mass_true'].SetLineColor(kRed)
    hists['res_right_vs_mass_true'] = TGraphErrors()
    hists['res_right_vs_mass_true'].SetName('res_right'+auxNameStr+'_vs_mass_true')
    hists['res_right_vs_mass_true'].SetMarkerSize(0.3)
    hists['res_right_vs_mass_true'].SetMarkerStyle(20)
    hists['res_right_vs_mass_true'].SetMarkerColor(kBlue)
    hists['res_right_vs_mass_true'].SetLineColor(kBlue)
    #-----------
    hists['mass_true']=TH1F('mass_true','m_{true}',200,0.,200.)

    reader = TMVA.Reader()

    for varName, var in sorted(variables.iteritems()):
        reader.AddVariable(varName,var[4])

    reader.BookMVA('BRT_HiggsMass',GetDir(_weightsDirName)+'/'+xmlFileName)

    resolution_at_125 = []

    min_calib_efficiency=0.95
    iPoint=0

    keys = []
    if not _useZtt:
        keys = mass_points
    else:
        for mass in mass_points:
            if mass < _Zmass:
                keys.append(mass)
        #Only calibrate on Higgs Samples, but test on Z
        if not(doCalibration):
            keys.append("Z")
        for mass in mass_points:
            if mass > _Zmass:
                keys.append(mass)

    fitmin = 60
    fitmax = 115
    if testCalibration:
        fitmin = 60
        fitmax = 170
    chi2=0.
    nInFitRange=0.
    for p, key in enumerate(keys):
        isZ = (key == "Z")
        if not isZ: mass = key
        if     isZ: mass = _Zmass
        print 'Mass:', mass
        nEvtsFail = 0
        nEvts = 0
        nEntries = trees[str(key)].GetEntriesFast()
        #nEntries=10000
        #Declare variables to calculate average mass of events
        sum_mass=0
        num_good_events=0
        for ientry in xrange(nEntries):
            if ientry%10000==0:
                print 'Event', ientry, 'out of', nEntries
            ## Get the next tree in the chain and verify
            if trees[str(key)].LoadTree(ientry) <  0: break
            ## Copy next entry into memory and verify
            if trees[str(key)].GetEntry(ientry) <= 0: continue
            ## Apply decay channel and selection cuts
            if not PassedCuts(False,True,trees[str(key)]): continue
            ## Load the variables used by the BRT
            for varName, var in sorted(variables.iteritems()):
                var[4][0] = getattr(trees[str(key)],varName)
                hists[varName].Fill(getattr(trees[str(key)],varName))
            ## Get the BRT output
            mass_BRT = reader.EvaluateMVA("BRT_HiggsMass")/1000.
            ## Use the calibration curve to correct the BRT output
            mass_BRT_corr = mass_BRT
            #Iterate variables to calculate average mass
            sum_mass+=trees[str(key)].ditau_m/1000.
            num_good_events+=1
            if testCalibration:
                nEvts += 1
                if calibration[0][1] <= mass_BRT <= calibration[-1][1]:
                    for c in range(len(calibration)):
                        if calibration[c][1] <= mass_BRT <= calibration[c+1][1]:
                            ss = (calibration[c+1][0]-calibration[c][0])/(calibration[c+1][1]-calibration[c][1])
                            mass_BRT_corr = calibration[c][0] + ss*(mass_BRT-calibration[c][1])
                else:
                    nEvtsFail += 1
                    continue
            mass_true = trees[str(key)].ditau_m/1000.
            hists['mass_true'].Fill(mass_true)
            hists['mass_BRT_'+str(key)].Fill(mass_BRT_corr)
            hists['mass_BRT_bias'].Fill(mass_BRT_corr-mass_true)
            if isZ: hists['mass_true_Z'].Fill(mass_true)
            hists['mass_true_all'].Fill(mass_true)
            hists['mass_BRT_all'].Fill(mass_BRT_corr)
            if fitmin<=trees[str(key)].ditau_m/1000.<=fitmax:
                chi2+=(mass_BRT_corr-trees[str(key)].ditau_m/1000.)**2
                nInFitRange+=1.
        if 125 <= mass <= 126:
            resolution_at_125 = [mass,rms/mean]
        average_mass=sum_mass/num_good_events
        if _useAverageMass:
            mass=average_mass
            print "Changed mass to", mass
        is_well_calibrated=True
        if testCalibration:
            if nEvts > 0:
                calib_efficiency = float((nEvts-nEvtsFail))/nEvts
            else:
                calib_efficiency = 0
            if isZ:
                print 'Calibration efficiency for Z sample = '+str(calib_efficiency*100)+'%'
            else:
                hists['calib_efficiency_vs_mass_true'].SetPoint(p,mass,calib_efficiency)
                hists['calib_efficiency_vs_mass_true'].SetPointError(p,0,0)
            if calib_efficiency<min_calib_efficiency:
                is_well_calibrated=False
        mean = hists['mass_BRT_'+str(key)].GetMean()
        mean_error = hists['mass_BRT_'+str(key)].GetMeanError()
        maximum = hists['mass_BRT_'+str(key)].GetXaxis().GetBinCenter(hists['mass_BRT_'+str(key)].GetMaximumBin())
        maximum_error = mean_error
        median = array.array('d',[0])
        yq = array.array('d',[0.5])
        hists['mass_BRT_'+str(key)].GetQuantiles(1,median,yq)
        median_error = mean_error
        rms = hists['mass_BRT_'+str(key)].GetRMS()
        rms_error = hists['mass_BRT_'+str(key)].GetRMSError()
        hists['mean_mass_BRT_vs_mass_true'].SetPoint(p,mass,mean)
        hists['mean_mass_BRT_vs_mass_true'].SetPointError(p,0,mean_error)
        hists['max_mass_BRT_vs_mass_true'].SetPoint(p,mass,maximum)
        hists['max_mass_BRT_vs_mass_true'].SetPointError(p,0,maximum_error)
        hists['med_mass_BRT_vs_mass_true'].SetPoint(p,mass,median[0])
        hists['med_mass_BRT_vs_mass_true'].SetPointError(p,0,median_error)
        hists['rms_mass_BRT_vs_mass_true'].SetPoint(p,mass,rms)
        hists['rms_mass_BRT_vs_mass_true'].SetPointError(p,0,rms_error)
        hists['res_mass_BRT_vs_mass_true'].SetPoint(p,mass,rms/mean)
        hists['res_mass_BRT_vs_mass_true'].SetPointError(p,0,math.sqrt(pow(rms_error/mean,2)+pow(rms/(mean*mean)*mean_error,2)))

        if doCalibration:
            calibration.append([mass,mean])
            

        #Fit Bifurcated Gaussian to peak
        Mass=RooRealVar("Mass","Invariant Mass",0.,200.)
        Sigma_left=RooRealVar("WidthL","Right Mass Width",5.,40.)
        Sigma_right=RooRealVar("WidthR","Right Mass Width",5.,40.)
        Mean=RooRealVar("Mean","Mean Invariant Mass",0.9*mean,1.1*mean)
        Mass.setRange("Subrange",0.75*mean,1.25*mean)
        
        FitData=RooDataHist("data","My Dataset",RooArgList(Mass),hists['mass_BRT_'+str(key)])
        
        BifurGauss=RooBifurGauss("bifurguass","Bifurcated Gaussian",Mass,Mean,Sigma_left,Sigma_right)
        myfit=BifurGauss.fitTo(FitData,RooFit.Save(),RooFit.Range("Subrange"))
        xframe=Mass.frame(40,200)
        FitData.plotOn(xframe,RooFit.Binning(20))
        BifurGauss.plotOn(xframe,RooFit.Range("Subrange"))
        xframe.SetTitle(' ')
        xframe.SetXTitle(mBRTleg+' [GeV]')
        xframe.SetYTitle('Entries')
        xframe.GetXaxis().SetTitleOffset(1.3)
        xframe.GetYaxis().SetTitleOffset(1.9)
        xframe.Draw()
        canvas.SaveAs(outdir+'/bifurgaussfit_mass_BRT_'+auxNameStr+str(key)+'.png')
        canvas.Clear()

        if is_well_calibrated:
            s_l=Sigma_left.getVal()
            s_r=Sigma_right.getVal()
            s_l_e=Sigma_left.getError()
            s_r_e=Sigma_right.getError()
            m=Mean.getVal()
            m_e=Mean.getError()
            
            hists['sigma_left_vs_mass_true'].SetPoint(iPoint,mass,s_l)
            hists['sigma_left_vs_mass_true'].SetPointError(iPoint,0,s_l_e)
            hists['sigma_right_vs_mass_true'].SetPoint(iPoint,mass,s_r)
            hists['sigma_right_vs_mass_true'].SetPointError(iPoint,0,s_r_e)
            hists['res_left_vs_mass_true'].SetPoint(iPoint,mass,s_l/m)
            hists['res_left_vs_mass_true'].SetPointError(iPoint,0,math.sqrt(pow(s_l_e/m,2)+pow(s_l/(m*m)*m_e,2)))
            hists['res_right_vs_mass_true'].SetPoint(iPoint,mass,s_r/m)
            hists['res_right_vs_mass_true'].SetPointError(iPoint,0,math.sqrt(pow(s_r_e/m,2)+pow(s_r/(m*m)*m_e,2)))
            iPoint+=1

    xmin = mass_points[0]-10
    xmax = mass_points[-1]+10
    ymin = mass_points[0]-10
    ymax = mass_points[-1]+10
    if doCalibration:
        pickle.dump(calibration,calibFile)

    ## Save all plots in the output file and as png 
    ## Save all plots in the output file and as png figures.
    #-----------
    if _useZtt and not(doCalibration):
        hists['mass_true_Z'].Draw()
        hists['mass_true_Z'].Write()
        canvas.SaveAs(outdir+'/'+hists['mass_true_Z'].GetName()+'.png')
        canvas.Clear()
    #-----------
    for key in keys:
        hists['mass_BRT_'+str(key)].Rebin(9)
        hists['mass_BRT_'+str(key)].Draw()
        hists['mass_BRT_'+str(key)].Write()
        canvas.SaveAs(outdir+'/'+hists['mass_BRT_'+str(key)].GetName()+'.png')
        canvas.Clear()
    #-----------
    hists['mass_true_all'].Draw()
    hists['mass_true_all'].Write()
    canvas.SaveAs(outdir+'/'+hists['mass_true_all'].GetName()+'.png')
    canvas.Clear()
    #-----------
    hists['mass_BRT_all'].Draw()
    hists['mass_BRT_all'].Write()
    canvas.SaveAs(outdir+'/'+hists['mass_BRT_all'].GetName()+'.png')
    canvas.Clear()
    #-----------
    hists['mass_BRT_bias'].Draw()
    hists['mass_BRT_bias'].Write()
    canvas.SaveAs(outdir+'/'+hists['mass_BRT_bias'].GetName()+'.png')
    canvas.Clear()
    #-----------
    fitmin = 60
    fitmax = 115
    if testCalibration:
        fitmin = 60
        fitmax = 170
    flin1 = TF1('flin1','[0]+[1]*x',fitmin-0.1,fitmax+0.1)
    flin1.SetLineColor(kRed)
    hists['mean_mass_BRT_vs_mass_true'].Fit(flin1,'RS','',fitmin-0.1,fitmax+0.1)
    intercept1 = flin1.GetParameter(0)
    slope1     = flin1.GetParameter(1)
    hists['mean_mass_BRT_vs_mass_true'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['mean_mass_BRT_vs_mass_true'].GetYaxis().SetTitle('<'+mBRTleg+'> [GeV]')
    hists['mean_mass_BRT_vs_mass_true'].GetXaxis().SetTitleOffset(1.3)
    hists['mean_mass_BRT_vs_mass_true'].GetYaxis().SetTitleOffset(1.9)
    hists['mean_mass_BRT_vs_mass_true'].GetXaxis().SetRangeUser(xmin,xmax)
    hists['mean_mass_BRT_vs_mass_true'].GetYaxis().SetRangeUser(ymin,ymax)
    hists['mean_mass_BRT_vs_mass_true'].Draw('AP')

    box.SetFillColor(kWhite)
    box.DrawBox((xmin+xmax)/2-5,ymin+(ymax-ymin)/16*4.5,xmax,ymin+(ymax-ymin)/16*7)
    latex.DrawLatex((xmin+xmax)/2,ymin+(ymax-ymin)/16*6,'Linear Fit in [%i,%i] GeV:'%(fitmin,fitmax))
    latex.DrawLatex((xmin+xmax)/2,ymin+(ymax-ymin)/16*5,mBRTleg+' = %0.3f#timesm_{true} + %0.1f'%(slope1,intercept1))
    latex.DrawLatex(xmin+(xmax-xmin)/16,ymin+(ymax-ymin)/16*14.5,channel_string)
    canvas.SaveAs(outdir+'/'+hists['mean_mass_BRT_vs_mass_true'].GetName()+'.png')
    canvas.Write('c_'+hists['mean_mass_BRT_vs_mass_true'].GetName())
    canvas.Clear()
    #-----------
    fitmin = 60
    fitmax = 115
    if testCalibration:
        fitmin = 60
        fitmax = 170
    flin2 = TF1('flin2','[0]+[1]*x',fitmin-0.1,fitmax+0.1)
    flin2.SetLineColor(kBlue)
    hists['max_mass_BRT_vs_mass_true'].Fit(flin2,'RS','',fitmin-0.1,fitmax+0.1)
    intercept2 = flin2.GetParameter(0)
    slope2     = flin2.GetParameter(1)
    hists['max_mass_BRT_vs_mass_true'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['max_mass_BRT_vs_mass_true'].GetYaxis().SetTitle('peak('+mBRTleg+') [GeV]')
    hists['max_mass_BRT_vs_mass_true'].GetXaxis().SetTitleOffset(1.3)
    hists['max_mass_BRT_vs_mass_true'].GetYaxis().SetTitleOffset(1.9)
    hists['max_mass_BRT_vs_mass_true'].GetXaxis().SetRangeUser(xmin,xmax)
    hists['max_mass_BRT_vs_mass_true'].GetYaxis().SetRangeUser(ymin,ymax)
    hists['max_mass_BRT_vs_mass_true'].Draw('AP')
    box.SetFillColor(kWhite)
    box.DrawBox((xmin+xmax)/2-5,ymin+(ymax-ymin)/16*4.5,xmax,ymin+(ymax-ymin)/16*7)
    latex.DrawLatex((xmin+xmax)/2,ymin+(ymax-ymin)/16*6,'Linear Fit in [%i,%i] GeV:'%(fitmin,fitmax))
    latex.DrawLatex((xmin+xmax)/2,ymin+(ymax-ymin)/16*5,mBRTleg+' = %0.3f#timesm_{true} + %0.1f'%(slope2,intercept2))
    latex.DrawLatex(xmin+(xmax-xmin)/16,ymin+(ymax-ymin)/16*14.5,channel_string)
    canvas.SaveAs(outdir+'/'+hists['max_mass_BRT_vs_mass_true'].GetName()+'.png')
    canvas.Write('c_'+hists['max_mass_BRT_vs_mass_true'].GetName())
    canvas.Clear()
    #-----------
    fitmin = 60
    fitmax = 115
    if testCalibration:
        fitmin = 60
        fitmax = 170
    flin3 = TF1('flin3','[0]+[1]*x',fitmin-0.1,fitmax+0.1)
    flin3.SetLineColor(kGreen)
    hists['med_mass_BRT_vs_mass_true'].Fit(flin3,'RS','',fitmin-0.1,fitmax+0.1)
    intercept3 = flin3.GetParameter(0)
    slope3     = flin3.GetParameter(1)
    hists['med_mass_BRT_vs_mass_true'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['med_mass_BRT_vs_mass_true'].GetYaxis().SetTitle('median('+mBRTleg+') [GeV]')
    hists['med_mass_BRT_vs_mass_true'].GetXaxis().SetTitleOffset(1.3)
    hists['med_mass_BRT_vs_mass_true'].GetYaxis().SetTitleOffset(1.9)
    hists['med_mass_BRT_vs_mass_true'].GetXaxis().SetRangeUser(xmin,xmax)
    hists['med_mass_BRT_vs_mass_true'].GetYaxis().SetRangeUser(ymin,ymax)
    hists['med_mass_BRT_vs_mass_true'].Draw('AP')
    
    box.SetFillColor(kWhite)
    box.DrawBox((xmin+xmax)/2-5,ymin+(ymax-ymin)/16*4.5,xmax,ymin+(ymax-ymin)/16*7)
    latex.DrawLatex((xmin+xmax)/2,ymin+(ymax-ymin)/16*6,'Linear Fit in [%i,%i] GeV:'%(fitmin,fitmax))
    latex.DrawLatex((xmin+xmax)/2,ymin+(ymax-ymin)/16*5,mBRTleg+' = %0.3f#timesm_{true} + %0.1f'%(slope3,intercept3))
    latex.DrawLatex(xmin+(xmax-xmin)/16,ymin+(ymax-ymin)/16*14.5,channel_string)
    canvas.SaveAs(outdir+'/'+hists['med_mass_BRT_vs_mass_true'].GetName()+'.png')
    canvas.Write('c_'+hists['med_mass_BRT_vs_mass_true'].GetName())
    canvas.Clear()
    #-----------
    ymin = 0.0
    ymax = 0.7
    hists['res_mass_BRT_vs_mass_true'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['res_mass_BRT_vs_mass_true'].GetYaxis().SetTitle('RMS('+mBRTleg+')/<'+mBRTleg+'>')
    hists['res_mass_BRT_vs_mass_true'].GetXaxis().SetTitleOffset(1.3)
    hists['res_mass_BRT_vs_mass_true'].GetYaxis().SetTitleOffset(2.0)
    hists['res_mass_BRT_vs_mass_true'].GetXaxis().SetRangeUser(xmin,xmax)
    hists['res_mass_BRT_vs_mass_true'].GetYaxis().SetRangeUser(ymin,ymax)
    hists['res_mass_BRT_vs_mass_true'].Draw('AP')
    line.SetLineStyle(kDashed)
    line.SetLineWidth(2)
    if len(resolution_at_125) == 2:
        line.DrawLine(resolution_at_125[0],0,resolution_at_125[0],resolution_at_125[1])
        line.DrawLine(xmin,resolution_at_125[1],resolution_at_125[0],resolution_at_125[1])
        latex.DrawLatex((xmin+xmax)/4,ymin+(ymax-ymin)/16*6,'#sigma_{m}/m at nominal %i GeV = %0.1f%%'%(resolution_at_125[0],resolution_at_125[1]*100))
    latex.DrawLatex(xmin+(xmax-xmin)/16,ymin+(ymax-ymin)/16*14.5,channel_string)
    canvas.SaveAs(outdir+'/'+hists['res_mass_BRT_vs_mass_true'].GetName()+'.png')
    canvas.Write('c_'+hists['res_mass_BRT_vs_mass_true'].GetName())
    canvas.Clear()
    #-----------
    ymin = 0.0
    ymax = 50.0
    hists['rms_mass_BRT_vs_mass_true'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['rms_mass_BRT_vs_mass_true'].GetYaxis().SetTitle('RMS('+mBRTleg+') [GeV]')
    hists['rms_mass_BRT_vs_mass_true'].GetXaxis().SetTitleOffset(1.3)
    hists['rms_mass_BRT_vs_mass_true'].GetYaxis().SetTitleOffset(1.8)
    hists['rms_mass_BRT_vs_mass_true'].GetXaxis().SetRangeUser(xmin,xmax)
    hists['rms_mass_BRT_vs_mass_true'].GetYaxis().SetRangeUser(ymin,ymax)
    hists['rms_mass_BRT_vs_mass_true'].Draw('AP')
    latex.DrawLatex(xmin+(xmax-xmin)/16,ymin+(ymax-ymin)/16*14.5,channel_string)
    canvas.SaveAs(outdir+'/'+hists['rms_mass_BRT_vs_mass_true'].GetName()+'.png')
    canvas.Write('c_'+hists['rms_mass_BRT_vs_mass_true'].GetName())
    canvas.Clear()
    #-----------
    if testCalibration:
        ymin = 0.0
        ymax = 1.2
        hists['calib_efficiency_vs_mass_true'].GetXaxis().SetTitle('m_{true} [GeV]')
        hists['calib_efficiency_vs_mass_true'].GetYaxis().SetTitle('calibration efficiency')
        hists['calib_efficiency_vs_mass_true'].GetXaxis().SetTitleOffset(1.3)
        hists['calib_efficiency_vs_mass_true'].GetYaxis().SetTitleOffset(1.8)
        hists['calib_efficiency_vs_mass_true'].GetXaxis().SetRangeUser(xmin,xmax)
        hists['calib_efficiency_vs_mass_true'].GetYaxis().SetRangeUser(ymin,ymax)
        hists['calib_efficiency_vs_mass_true'].Draw('AP')
        latex.DrawLatex(xmin+(xmax-xmin)/16,ymin+(ymax-ymin)/16*14.5,channel_string)
        canvas.SaveAs(outdir+'/'+hists['calib_efficiency_vs_mass_true'].GetName()+'.png')
        canvas.Write('c_'+hists['calib_efficiency_vs_mass_true'].GetName())
        canvas.Clear()
    #-----------
    for varName, var in sorted(variables.iteritems()):
        hists[varName].Write()
    #-----------
    #if doCalibration or testCalibration:
    outFile.Close()

    for key in keys:
        files[str(key)].Close()

    if doCalibration:
        calibFile.close()
    elif testCalibration:
        calibFile.close()

    #----------
    ymin = 0.0
    ymax = 50.0
    hists['sigma_left_vs_mass_true'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['sigma_left_vs_mass_true'].GetYaxis().SetTitle('RMS('+mBRTleg+') [GeV]')
    hists['sigma_left_vs_mass_true'].GetXaxis().SetTitleOffset(1.3)
    hists['sigma_left_vs_mass_true'].GetYaxis().SetTitleOffset(1.8)
    hists['sigma_left_vs_mass_true'].GetXaxis().SetRangeUser(xmin,xmax)
    hists['sigma_left_vs_mass_true'].GetYaxis().SetRangeUser(ymin,ymax)
    hists['sigma_left_vs_mass_true'].Draw('AP')
    #--------
    hists['sigma_right_vs_mass_true'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['sigma_right_vs_mass_true'].GetYaxis().SetTitle('RMS('+mBRTleg+') [GeV]')
    hists['sigma_right_vs_mass_true'].GetXaxis().SetTitleOffset(1.3)
    hists['sigma_right_vs_mass_true'].GetYaxis().SetTitleOffset(1.8)
    hists['sigma_right_vs_mass_true'].GetXaxis().SetRangeUser(xmin,xmax)
    hists['sigma_right_vs_mass_true'].GetYaxis().SetRangeUser(ymin,ymax)
    hists['sigma_right_vs_mass_true'].Draw('P same')
    legend=TLegend(0.7,0.7,0.9,0.9)
    legend.AddEntry(hists['sigma_left_vs_mass_true'],'Left RMS','p')
    legend.AddEntry(hists['sigma_right_vs_mass_true'],'Right RMS','p')
    legend.Draw('same')
    latex.DrawLatex(xmin+(xmax-xmin)/16,ymin+(ymax-ymin)/16*14.5,channel_string)
    canvas.SaveAs(outdir+'/sigma_vs_mass_true'+auxNameStr+'.png')
    canvas.Clear()
        
    #----------
    ymax = 0.7
    xmin=hists['res_left_vs_mass_true'].GetX()[0]
    print xmin
    hists['res_left_vs_mass_true'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['res_left_vs_mass_true'].GetYaxis().SetTitle('RMS('+mBRTleg+')/<'+mBRTleg+'>')
    hists['res_left_vs_mass_true'].GetXaxis().SetTitleOffset(1.3)
    hists['res_left_vs_mass_true'].GetYaxis().SetTitleOffset(1.8)
    hists['res_left_vs_mass_true'].GetXaxis().SetRangeUser(xmin,xmax)
    hists['res_left_vs_mass_true'].GetYaxis().SetRangeUser(ymin,ymax)
    hists['res_left_vs_mass_true'].Draw('AP')
    #--------
    hists['res_right_vs_mass_true'].GetXaxis().SetTitle('m_{true} [GeV]')
    hists['res_right_vs_mass_true'].GetYaxis().SetTitle('RMS('+mBRTleg+')/<'+mBRTleg+'>')
    hists['res_right_vs_mass_true'].GetXaxis().SetTitleOffset(1.3)
    hists['res_right_vs_mass_true'].GetYaxis().SetTitleOffset(1.8)
    hists['res_right_vs_mass_true'].GetXaxis().SetRangeUser(xmin,xmax)
    hists['res_right_vs_mass_true'].GetYaxis().SetRangeUser(ymin,ymax)
    hists['res_right_vs_mass_true'].Draw('P same')
    legend=TLegend(0.7,0.7,0.9,0.9)
    legend.AddEntry(hists['res_left_vs_mass_true'],'Left Resolution','p')
    legend.AddEntry(hists['res_right_vs_mass_true'],'Right Resolution','p')
    legend.Draw('same')
    latex.DrawLatex(xmin+(xmax-xmin)/16,ymin+(ymax-ymin)/16*14.5,channel_string)
    canvas.SaveAs(outdir+'/res_vs_mass_true'+auxNameStr+'.png')
    canvas.Clear()
    
    return slope1


#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
## Auxiliary functions used to apply cuts in two different ways:                                   ##
## 1) using a string in a TCut class --> this is the case used in training                         ##
## 2) applying the cuts using the trees variables --> this is the case used in testing             ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################


def PassedCuts(asString,YesOrNo,tree=TTree()):
    return Cuts(asString,YesOrNo,tree)

def Cuts(asString,YesOrNo,tree=TTree()):
    if not OR(asString,YesOrNo):
        print 'In Cuts(): Error: one of the input variables "asString" or "YesOrNo" must be = True.'
        sys.exit()
    if not XOR(asString,YesOrNo):
        print 'In Cuts(): Error: input variables "asString" and "YesOrNo" can not be both = True.'
        sys.exit()
    if _useAnalysisPresel:
        met_min=20000
        lep1_pt_min=35000
        lep2_pt_min=25000
        lep1_eta_max=2.5
        lep2_eta_max=2.5
        dR_lep_lep_max=2.6
        dR_lep_lep_min=0.8
    else:
        met_min = 20000.
        lep1_pt_min = 15000.
        lep2_pt_min = 15000.
        lep1_eta_max = 2.4
        lep2_eta_max = 2.4
        dR_lep_lep_max=math.sqrt(math.pi**2+(lep1_eta_max+lep2_eta_max)**2)
        dR_lep_lep_min=0.
    if YesOrNo:
        passed = True
        if _doLepLep and tree.leplep == 0: passed = False
        if _doLepHad and tree.lephad == 0: passed = False
        if _doHadHad and tree.hadhad == 0: passed = False
        if _applyPreselCuts:
            if _useAnalysisVariables:
                if passed and not (abs(tree.lep1_eta) < lep1_eta_max): passed = False
                if passed and not (abs(tree.lep2_eta) < lep2_eta_max): passed = False
                if passed and not (tree.lep1_pt > lep1_pt_min): passed = False
                if passed and not (tree.lep2_pt > lep2_pt_min): passed = False
                if passed and not (tree.met_et > met_min): passed = False
                if passed and not (dR_lep_lep_min<tree.dR_lep_lep<dR_lep_lep_max): passed = False
            if _useSmearedAnalysisVariables:
                if passed and not (abs(tree.lep1_eta_sm) < lep1_eta_max): passed = False
                if passed and not (abs(tree.lep2_eta_sm) < lep2_eta_max): passed = False
                if passed and not (tree.lep1_pt_sm > lep1_pt_min): passed = False
                if passed and not (tree.lep2_pt_sm > lep2_pt_min): passed = False
                if passed and not (tree.met_et_sm > met_min): passed = False
                if passed and not (dR_lep_lep_min<tree.dR_lep_lep_sm<dR_lep_lep_max): passed = False
    if asString:
        decay_channel_cuts = ''
        if _doLepLep: decay_channel_cuts = 'leplep==1'
        if _doLepHad: decay_channel_cuts = 'lephad==1'
        if _doHadHad: decay_channel_cuts = 'hadhad==1'
        kinematic_cuts = ''
        if _applyPreselCuts:
            if _useAnalysisVariables:
                kinematic_cuts  = '    abs(lep1_eta) < '+str(lep1_eta_max)
                kinematic_cuts += ' && abs(lep2_eta) < '+str(lep2_eta_max)
                kinematic_cuts += ' && lep1_pt > '+str(lep1_pt_min)
                kinematic_cuts += ' && lep2_pt > '+str(lep2_pt_min)
                kinematic_cuts += ' && met_et > '+str(met_min)
            if _useSmearedAnalysisVariables:
                kinematic_cuts  = '    abs(lep1_eta_sm) < '+str(lep1_eta_max)
                kinematic_cuts += ' && abs(lep2_eta_sm) < '+str(lep2_eta_max)
                kinematic_cuts += ' && lep1_pt_sm > '+str(lep1_pt_min)
                kinematic_cuts += ' && lep2_pt_sm > '+str(lep2_pt_min)
                kinematic_cuts += ' && met_et_sm > '+str(met_min)
        cuts = decay_channel_cuts+' && '+kinematic_cuts
        while cuts[0]  in [' ','&','|']: cuts = cuts[1:]
        while cuts[-1] in [' ','&','|']: cuts = cuts[:-2]
        return cuts


#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
## Function used to return the (sub)directory where the input files are, where to write            ##
## the weights file, and where to put the plots and root files.                                    ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################


def GetDir(basedir):
    directory = basedir
    if 'trees' in basedir:
        if _useBosonVariables:           directory += '/with_boson_variables/'
        if _useTauVariables:             directory += '/with_tau_variables/'
        if _useTauDecayVariables:        directory += '/with_tau_decay_variables/'
        if _useAnalysisVariables:        directory += '/with_analysis_variables/'
        if _useSmearedAnalysisVariables: directory += '/with_smeared_analysis_variables/'
    else:
        if not(_target=='ditau_m'): directory+='/'+_target+'/'
        if _useZtt:            directory += '/with_Z/'
        if not(_useAverageMass):  directory += '/nominal_mass/'
        if _applyPreselCuts and _useAnalysisPresel:   directory += '/applying_analysis_presel_cuts/'
        elif _applyPreselCuts: directory += '/applying_preselection_cuts/'
        else:                 directory += '/without_preselection_cuts/'
        if _useBosonVariables: directory += '/with_boson_variables/'
        if _useTauVariables:   directory += '/with_tau_variables/'
        if _useTauDecayVariables:
            if _doLepLep:      directory += '/with_tau_decay_variables/leplep/'
            if _doLepHad:      directory += '/with_tau_decay_variables/lephad/'
            if _doHadHad:      directory += '/with_tau_decay_variables/hadhad/'
        if _useAnalysisVariables:
            if _doLepLep:      directory += '/with_analysis_variables/leplep/'
            if _doLepHad:      directory += '/with_analysis_variables/lephad/'
            if _doHadHad:      directory += '/with_analysis_variables/hadhad/'
        if _useSmearedAnalysisVariables:
            if _doLepLep:      directory += '/with_smeared_analysis_variables/leplep/'
            if _doLepHad:      directory += '/with_smeared_analysis_variables/lephad/'
            if _doHadHad:      directory += '/with_smeared_analysis_variables/hadhad/'
    return directory.replace('//','/')



##----------------------------------------------------------------------------------------


def inputHttSample(mass):
    return _HttSamplesDirName+'/Hmass%s.root'%(mass)
def HttSample(mass):
    return GetDir(_treesDirName)+'/H/ggHtautau_m%s.root'%(mass)

def inputZttSample():
    return _ZttSamplesDirName+'/Z.root'
def ZttSample():
    return GetDir(_treesDirName)+'/Z/Ztautau.root'


#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
## Function to calculate the Met_Phi_Centrality variable, as well as some basic phi operations.    ## 
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################

def getMetPhiCentrality(Tau1_phi,Tau2_phi,MET_phi):
    A=math.sin(MET_phi-Tau1_phi)/math.sin(Tau2_phi-Tau1_phi)
    B=math.sin(Tau2_phi-MET_phi)/math.sin(Tau2_phi-Tau1_phi)
    return (A+B)/math.sqrt(A**2+B**2)

def deltaPhi(phi1,phi2):
    dPhi=abs(phi1-phi2)
    while dPhi>math.pi:
        dPhi=dPhi-2*math.pi
    return abs(dPhi)

def addPhi(phi1,phi2):
    sumPhi=phi1+phi2
    while sumPhi>math.pi:
        sumPhi-=2*math.pi
    while sumPhi<-math.pi:
        sumPhi+=2*math.pi
    return sumPhi

#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
## Auxiliary function to read user input when program starts to run. The idea was to use           ##
## this function to confirm that the user wants to run with the given settings. But I also         ##
## want that if the user doesn't give an answer after some time, then code should run              ##
## anyway, because this means that the user sent the jobs and went to sleep. I couldn't            ##
## make this last part to work, so this function is not used so far.                               ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################

# NOT USED
def input_with_timeout(seconds,text):
    answer = None
    def timeout():
        subprocess.Popen(['/bin/bash','-c','echo timeout'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        answer = 'timeout'
    def countdown():
        remaining_time = seconds
        while remaining_time > 0:
            sys.stdout.write('\rRemaining time %i seconds ...'%(remaining_time))
            sys.stdout.flush()
            sleep(1)
            remaining_time -= 1
        sys.stdout.write('\n')
    timeout_timer = Timer(seconds,timeout)
    timeout_timer.start()
    try:
        print text,
        answer = sys.stdin.readline()
    except KeyboardInterrupt:
        answer = None
    timeout_timer.cancel()
    return answer


#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
##                                Here is where the code runs                                      ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################


slopes_bef_calib = {}
slopes_aft_calib = {}

if _fillTrainTestTrees:
    print 'Will save train and test trees in '+GetDir(_treesDirName)
    for mass in _mass_points_all:
        iFileName = inputHttSample(mass)
        oFileName = HttSample(mass)
        nTrainEvtsMax = 50000
        FillTrainTestTrees(iFileName.replace('//','/'),oFileName.replace('//','/'),nTrainEvtsMax)
    if _useZtt:
        iFileName = inputZttSample()
        oFileName = ZttSample()
        nTrainEvtsMax = 500000
        FillTrainTestTrees(iFileName.replace('//','/'),oFileName.replace('//','/'),nTrainEvtsMax)

for nem in [20]: #[10, 20, 40, 60, 80]:
    for nc in [20]: #-1 is not supported yet, the system sets it to 200
        for md in [100]:
            for nt in [40]: #[10, 20, 40, 60, 80, 100, 150, 200, 250]:
                factoryName = 'TMVARegression'
                methodName  = 'BRT_HiggsMass'
                trainParams = {}
                trainParams['AdaBoostBeta'] = 0.2
                trainParams['nCuts']        = nc
                trainParams['NTrees']       = nt
                trainParams['nEventsMin']   = nem
                trainParams['MaxDepth']     = md
                if _doTraining:
                    Training(factoryName,methodName,trainParams)
                if _doTesting:
                    xmlFileName = factoryName+'_'+methodName
                    training_parameters = ''
                    for key, value in trainParams.iteritems():
                        xmlFileName += '_'+str(key)+str(value)
                        training_parameters += '   '+str(key)+'='+str(value)
                    xmlFileName += '.weights.xml'
                    if _testWOCalibration:
                        slope = Testing(xmlFileName,True,False,False)
                        slopes_bef_calib[training_parameters] = slope
                    if _doCalibration:
                        slope = Testing(xmlFileName,False,True,False)
                        slopes_bef_calib[training_parameters] = slope
                    if _testCalibration:
                        slope = Testing(xmlFileName,False,False,True)
                        slopes_aft_calib[training_parameters] = slope

## Print information.
if _doTesting:
    if OR(_testWOCalibration,_doCalibration):
        print '---------------------------------------------------------------------------------------------------------------------------------------'
        print '|                             Training parameters                                | Slope of m_BRT vs m_true linear fit (before calib) |'
        print '---------------------------------------------------------------------------------------------------------------------------------------'
        for training_parameters, slope in slopes_bef_calib.iteritems():
            print '|%s       |           %s                           |'%(training_parameters,slope)
        print '---------------------------------------------------------------------------------------------------------------------------------------'
    if _testCalibration:
        print '---------------------------------------------------------------------------------------------------------------------------------------'
        print '|                             Training parameters                                | Slope of m_BRT vs m_true linear fit (after calib)  |'
        print '---------------------------------------------------------------------------------------------------------------------------------------'
        for training_parameters, slope in slopes_aft_calib.iteritems():
            print '|%s       |           %s                           |'%(training_parameters,slope)
        print '---------------------------------------------------------------------------------------------------------------------------------------'


#####################################################################################################
##-------------------------------------------------------------------------------------------------##
##                                                                                                 ##
##                                  ######   ##    ##   #####                                      ##
##                                  ##       ###   ##   ##   ##                                    ##
##                                  ######   ## ## ##   ##   ##                                    ##
##                                  ##       ##   ###   ##   ##                                    ##
##                                  ######   ##    ##   #####                                      ##
##                                                                                                 ##
##-------------------------------------------------------------------------------------------------##
#####################################################################################################
