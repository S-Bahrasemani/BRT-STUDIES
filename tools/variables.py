import math
from .mass import collinearmass
from . import tau
import ROOT

def getMetPhiCentrality(Tau1_phi,Tau2_phi,MET_phi):
    A=math.sin(MET_phi-Tau1_phi)/math.sin(Tau2_phi-Tau1_phi)
    B=math.sin(Tau2_phi-MET_phi)/math.sin(Tau2_phi-Tau1_phi)
    return (A+B)/math.sqrt(A**2+B**2)

def deltaPhi(phi1,phi2):
    dPhi=abs(phi1-phi2)
    while dPhi>math.pi:
        dPhi=dPhi-2*math.pi
    return abs(dPhi)

def evalVariable(var_name,Tau1,Tau2,met_px,met_py,Energy_scale=''):
    #When adding new variables, they must be added to this function
    energyVarList=['lep1_pt','lep2_pt','met_et','transverse_mass_lep1_lep2','transverse_mass_lep1_met','transverse_mass_lep2_met','ptsum_lep1_lep2_met','ptsum_lep1_lep2','collinear_mass','visible_mass','averaged_mass','resonance_pt','proj_mass_lep1','proj_mass_lep2']
    if 'smirnov_par_perp' in var_name:
        var_name=var_name.replace('_smirnov_par_perp','')
    if Energy_scale=='' or not(var_name in energyVarList):
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
	if var_name=='met_phi':
	    return math.atan2(met_py,met_px)
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
        if var_name=='pt_ratio_lep2_lep1':
            return Tau2.pt/Tau1.pt
        if var_name=='pt_ratio_met_ditau':
            return math.sqrt(met_px**2+met_py**2)/math.sqrt((Tau1.px+Tau2.px)**2+(Tau1.py+Tau2.py)**2)
        if var_name=='visible_mass':
            return (Tau1.tlv+Tau2.tlv).M()
        if var_name=='averaged_mass':
            MET_tlv=ROOT.TLorentzVector()
            MET_tlv.SetPtEtaPhiM(math.sqrt(met_px**2+met_py**2),(Tau1.eta+Tau2.eta)/2.,math.atan2(met_py,met_px),0.)
            return(Tau1.tlv+Tau2.tlv+MET_tlv).M()
        if var_name=='mass_fitted':
            return 1.936836*evalVariable('visible_mass',Tau1,Tau2,met_px,met_py)/1000.-48.6240
	if var_name=='resonance_pt':
	    return evalVariable('pttot_lep1_lep2_met',Tau1,Tau2,met_px,met_py)*evalVariable('ptsum_lep1_lep2_met',Tau1,Tau2,met_px,met_py)
        if var_name=='proj_mass_lep1':
            MET_tlv=ROOT.TLorentzVector()
            MET_tlv.SetPtEtaPhiM(math.sqrt(met_px**2+met_py**2),(Tau1.eta),math.atan2(met_py,met_px),0.)
            return (Tau1.tlv+Tau2.tlv+MET_tlv).M()
        if var_name=='proj_mass_lep2':
            MET_tlv=ROOT.TLorentzVector()
            MET_tlv.SetPtEtaPhiM(math.sqrt(met_px**2+met_py**2),(Tau2.eta),math.atan2(met_py,met_px),0.)
            return (Tau1.tlv+Tau2.tlv+MET_tlv).M()
        if var_name=='proj_mass_lep1_by_visible_mass':
            return evalVariable('proj_mass_lep1',Tau1,Tau2,met_px,met_py,'')/evalVariable('visible_mass',Tau1,Tau2,met_px,met_py,'')
        if var_name=='proj_mass_lep2_by_visible_mass':
            return evalVariable('proj_mass_lep2',Tau1,Tau2,met_px,met_py)/evalVariable('visible_mass',Tau1,Tau2,met_px,met_py,'')
        if var_name=='pttot_lep1_lep2_met*ptsum_lep1_lep2_met':
            return evalVariable('pttot_lep1_lep2_met',Tau1,Tau2,met_px,met_py)/evalVariable('ptsum_lep1_lep2_met',Tau1,Tau2,met_px,met_py,'')
        if var_name=='pttot_lep1_lep2*ptsum_lep1_lep2':
            return evalVariable('pttot_lep1_lep2',Tau1,Tau2,met_px,met_py)/evalVariable('ptsum_lep1_lep2',Tau1,Tau2,met_px,met_py,'')
        if var_name=='ptdiff_lep1_lep2*ptsum_lep1_lep2':
            return evalVariable('ptdiff_lep1_lep2',Tau1,Tau2,met_px,met_py)/evalVariable('ptsum_lep1_lep2',Tau1,Tau2,met_px,met_py,'')
        if var_name=='abs(lep1_eta)':
            return abs(Tau1.eta)
        if var_name=='lep2_eta*lep1_eta/abs(lep1_eta)':
            if Tau1.eta<0:
                return -Tau2.eta
            else:
                return Tau2.eta
        if var_name=='lep1_eta*lep2_eta':
            return Tau1.eta*Tau2.eta
        if var_name=='(lep2_phi-lep1_phi)%(2*3.141592)':
            return (Tau2.phi-Tau1.phi)%(2*3.141592)
        if var_name=='(met_phi-lep1_phi)%(2*3.141592)':
            return (math.atan2(met_py,met_px)-Tau1.phi)%(2*3.141592)
        print "Error, unknown variable name. Variable must be added to evalVariable function."
        return 0
    else:
        return evalVariable(var_name,Tau1,Tau2,met_px,met_py)/evalVariable(Energy_scale,Tau1,Tau2,met_px,met_py)

def addPhi(phi1,phi2):
    sumPhi=phi1+phi2
    while sumPhi>math.pi:
        sumPhi-=2*math.pi
    while sumPhi<-math.pi:
        sumPhi+=2*math.pi
    return sumPhi
