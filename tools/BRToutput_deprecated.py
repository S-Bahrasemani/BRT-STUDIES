from ROOT import TBranch, TMVA, TCanvas, TH1F 
import math
import array
import pickle
from .mass import collinearmass

withRelativeVariables=False
withRecoTraining=True

Energy_scale='ptsum_lep1_lep2_met'
energyvariablelist=['lep1_pt','lep2_pt','met_et','transverse_mass_lep1_lep2','transverse_mass_lep1_met','transverse_mass_lep2_met','ptsum_lep1_lep2_met','ptsum_lep1_lep2']

weightFilePath='/cluster/data04/mquennev/higgs/weights/alpha_vis/applying_preselection_cuts/with_full_sim_variables/hadhad/TMVARegression_BRT_HiggsMass_nEventsMin20_AdaBoostBeta0.2_nCuts20_MaxDepth100_NTrees40.weights.xml'

variables = {}
variables['lep1_pt']                   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
variables['lep1_eta']                  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
variables['lep2_pt']                   = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
variables['lep2_eta']                  = [''   ,'F',-10.00,10.00,  array.array('f',[0]),TBranch(), -10.00,10.00,  200]
variables['met_et']                    = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,200000., 40]
variables['transverse_mass_lep1_lep2'] = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 0]
variables['transverse_mass_lep1_met']  = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 50]
variables['transverse_mass_lep2_met']  = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 50]
variables['dphi_lep1_met']             = ['rad','F',  0.00, 3.15,  array.array('f',[0]),TBranch(),   0.00, 3.15,   30]
variables['dphi_lep2_met']             = ['rad','F',  0.00, 3.15,  array.array('f',[0]),TBranch(),   0.00, 3.15,   30]
variables['dphi_lep_lep']              = ['rad','F',  0.00, 3.15,  array.array('f',[0]),TBranch(),   0.00, 3.15,   30]
variables['deta_lep_lep']              = ['',   'F',  0.00,20.00,  array.array('f',[0]),TBranch(),   0.00,20.00,  200]
variables['dR_lep_lep']                = ['',   'F',  0.00,25.00,  array.array('f',[0]),TBranch(),   0.00,20.00,  200]
variables['ptsum_lep1_lep2']           = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,250000., 25]
variables['ptsum_lep1_lep2_met']       = ['MeV','F',  0.00,9999999,array.array('f',[0]),TBranch(),   0.00,300000., 30]
variables['pttot_lep1_lep2_met']       = ['',   'F',  0.00, 2.00,  array.array('f',[0]),TBranch(),   0.00, 1.10,   22]
variables['pttot_lep1_lep2']           = ['',   'F',  0.00, 2.00,  array.array('f',[0]),TBranch(),   0.00, 1.10,   22]
variables['ptdiff_lep1_lep2']          = ['',   'F',  0.00, 2.00,  array.array('f',[0]),TBranch(),   0.00, 1.10,   22]

reader=TMVA.Reader()
hists={}

for varName, var in sorted(variables.iteritems()):
    if withRecoTraining and withRelativeVariables and varName==Energy_scale:
        continue
    reader.AddVariable(varName,var[4])

if withRecoTraining:
    a=array.array('f',[0])
    reader.AddSpectator("truth_mass",a)
    reader.AddSpectator("visible_mass",a)
    reader.AddSpectator("mmc0_resonance_m",a)
    reader.AddSpectator("mmc1_resonance_m",a)
    reader.AddSpectator("mass_collinear_tau1_tau2",a)

reader.BookMVA('BRT_HiggsMass',weightFilePath)

def mass_BRT(Tau1,Tau2,met_px,met_py):
    UpdateVariables(Tau1,Tau2,met_px,met_py)
    return reader.EvaluateMVA("BRT_HiggsMass")

def deltaPhi(phi1,phi2):
    dphi=abs(phi1-phi2)
    while dphi>math.pi:
        dphi=dphi-2*math.pi
    return abs(dphi)

def getMetPhiCentrality(Tau1_phi,Tau2_phi,MET_phi):
    A=math.sin(MET_phi-Tau1_phi)/math.sin(Tau2_phi-Tau1_phi)
    B=math.sin(Tau2_phi-MET_phi)/math.sin(Tau2_phi-Tau1_phi)
    return (A+B)/math.sqrt(A**2+B**2)

def UpdateVariables(Tau1,Tau2,met_px,met_py):
    for varName, var in sorted(variables.iteritems()):
        if varName in energyvariablelist and withRelativeVariables:
            var[4][0] = evalVariable(varName,Tau1,Tau2,met_px,met_py)/Energy_scale
        else:
            var[4][0] = evalVariable(varName,Tau1,Tau2,met_px,met_py)

def evalVariable(var_name,Tau1,Tau2,met_px,met_py):
    #When adding new variables, they must be added to this functio
    if var_name=='lep1_pt':
        return Tau1.pt
    if var_name=='lep1_eta':
        return Tau1.eta
    if var_name=='lep2_pt':
        return Tau2.pt
    if var_name=='lep2_eta':
        return Tau2.eta
    if var_name=='met_et':
        return math.sqrt(met_px**2+met_py**2)
    if var_name=='transverse_mass_lep1_lep2':
        return math.sqrt(2*Tau1.pt*Tau2.pt*(1.-math.cos(deltaPhi(Tau1.phi,Tau2.phi))))
    if var_name=='transverse_mass_lep1_met':
        return math.sqrt(2*Tau1.pt*math.sqrt(met_px**2+met_py**2)*(1.-math.cos(deltaPhi(Tau1.phi,math.atan2(met_py,met_px)))))   
    if var_name=='transverse_mass_lep2_met':
        return math.sqrt(2*Tau2.pt*math.sqrt(met_px**2+met_py**2)*(1.-math.cos(deltaPhi(Tau2.phi,math.atan2(met_py,met_px)))))
    if var_name=='dphi_lep1_met':
        return deltaPhi(Tau1.phi,math.atan2(met_py,met_px))
    if var_name=='dphi_lep2_met':
        return deltaPhi(Tau2.phi,math.atan2(met_py,met_px))
    if var_name=='dphi_lep_lep':
        return deltaPhi(Tau1.phi,Tau2.phi)
    if var_name=='deta_lep_lep':
        return abs(Tau1.eta-Tau2.eta)
    if var_name=='dR_lep_lep':
        return math.sqrt((Tau1.eta-Tau2.eta)**2+deltaPhi(Tau1.phi,Tau2.phi)**2)
    if var_name=='ptsum_lep1_lep2_met':
        return Tau1.pt+Tau2.pt+math.sqrt(met_px**2+met_py**2)
    if var_name=='ptsum_lep1_lep2':
        return Tau1.pt+Tau2.pt
    if var_name=='pttot_lep1_lep2_met':
        return math.sqrt((Tau1.px+Tau2.px+met_px)**2+(Tau1.py+Tau2.py+met_py)**2)/(Tau1.pt+Tau2.pt+math.sqrt(met_px**2+met_py**2))
    if var_name=='pttot_lep1_lep2':
        return math.sqrt((Tau1.px+Tau2.px)**2+(Tau1.py+Tau2.py)**2)/(Tau1.pt+Tau2.pt)
    if var_name=='ptdiff_lep1_lep2':
        return math.sqrt((Tau1.px-Tau2.px)**2+(Tau1.py-Tau2.py)**2)/(Tau1.pt+Tau2.pt)
    if var_name=='met_phi_centrality':
        return getMetPhiCentrality(Tau1.phi,Tau2.phi,math.atan2(met_py,met_px))
    if var_name=='collinear_mass':
        return collinearmass.mass(Tau1,Tau2,met_px,met_py)[1]
    if var_name=='visible_mass':
        return (Tau1.tlv+Tau2.tlv).M()
    print "Error, unknown variable name. Variable must be added to evalVariable function."
    return 0
