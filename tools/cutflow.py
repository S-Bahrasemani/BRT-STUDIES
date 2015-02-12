import ROOT
import math
# All basic cut definitions are here

def NOT(a):
    return not(a)

def OR(a,b):
    return (a or b)

def AND(a,b):
    return (a and b)

def list_AND(cutlist):
    if len(cutlist)<2:
        print 'Error, "AND"ing less than two cuts'
        return True
    elif len(cutlist)==2:
        return AND(cutlist[0],cutlist[1])
    else:
        return AND(cutlist[0],list_AND(cutlist[1:]))

def passedCut(event,cutstring):
    if not(cutstring in ['preselection','boosted','vbf']):
        print 'Cutstring not found'
        return False

    TAU1_MEDIUM = event.tau1_JetBDTSigMedium==1
    TAU2_MEDIUM = event.tau2_JetBDTSigMedium==1
    TAU1_TIGHT = event.tau1_JetBDTSigTight==1
    TAU2_TIGHT = event.tau2_JetBDTSigTight==1

    ID_MEDIUM = AND(TAU1_MEDIUM,TAU2_MEDIUM)
    ID_TIGHT = AND(TAU1_TIGHT,TAU2_TIGHT)
    ID_MEDIUM_TIGHT = OR(AND(TAU1_MEDIUM,TAU2_TIGHT),AND(TAU1_TIGHT,TAU2_MEDIUM))
    # ID cuts for control region where both taus are medium but not tight
    ID_MEDIUM_NOT_TIGHT = AND(AND(TAU1_MEDIUM,NOT(TAU1_TIGHT)),AND(TAU2_MEDIUM,NOT(TAU2_TIGHT)))

    TAU_SAME_VERTEX = event.tau_same_vertex

    LEAD_TAU_35 = event.tau1_pt > 35000
    SUBLEAD_TAU_25 = event.tau2_pt > 25000

    LEAD_JET_50 = event.jet1_pt > 50000
    SUBLEAD_JET_30 = event.jet2_pt > 30000
    AT_LEAST_1JET = event.jet1_pt > 30000

    CUTS_2J = AND(LEAD_JET_50,SUBLEAD_JET_30)
    CUTS_1J = AND(LEAD_JET_50,NOT(SUBLEAD_JET_30))
    CUTS_0J = NOT(LEAD_JET_50)

    ETA_TAUS=AND(abs(event.tau1_eta)<2.4,abs(event.tau2_eta)<2.4)

    MET = event.MET_et > 20000
    DR_TAUS = 0.8 < event.dR_tau1_tau2 < 2.4
    DETA_TAUS = event.dEta_tau1_tau2 < 1.5
    RESONANCE_PT = event.resonance_pt > 100000

    DETA_JETS = event.dEta_jets>2.0

    MMC_MASS = event.mmc0_resonance_m>0.

    MET_CENTRALITY = event.MET_bisecting or (event.dPhi_min_tau_MET < math.pi/4.)

    # common preselection cuts
    PRESELECTION = (
        list_AND([LEAD_TAU_35,SUBLEAD_TAU_25,ID_MEDIUM_TIGHT,MET,MMC_MASS,DR_TAUS,TAU_SAME_VERTEX,MET_CENTRALITY])
    )

    # VBF category cuts
    CUTS_VBF = (
        list_AND([PRESELECTION,CUTS_2J,DETA_TAUS,DETA_JETS])
    )

    # Boosted category cuts
    CUTS_BOOSTED = (
        list_AND([PRESELECTION,NOT(CUTS_VBF),RESONANCE_PT,DETA_TAUS])
    )
    if cutstring=='preselection':
        return PRESELECTION
    elif cutstring=='vbf':
        return CUTS_VBF
    elif cutstring=='boosted':
        return CUTS_BOOSTED

def passedCut_truth(event):
    LEAD_TAU_35 = event.lep1_pt > 35000
    SUBLEAD_TAU_25 = event.lep2_pt > 25000

    MET = event.met_et > 20000
    DR_TAUS = 0.8 < event.dR_lep_lep < 2.4

    MET_CENTRALITY = event.is_met_phi_central or (min(event.dphi_lep1_met,event.dphi_lep2_met) < math.pi/4.)

    # common preselection cuts
    PRESELECTION = (
        list_AND([LEAD_TAU_35,SUBLEAD_TAU_25,MET,DR_TAUS,MET_CENTRALITY])
    )
    return PRESELECTION
