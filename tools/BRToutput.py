import ROOT
import array
import variables

class BRTMass:
    def __init__(self,weightFilePath,useCorrection=False,Corrected_variable='',useRelativeVariables=False,Energy_scale=''):
	self.useCorrection=useCorrection
	self.correctionFactor=Corrected_variable
        self.withRelativeVariables=useRelativeVariables
	self.Energy_scale=Energy_scale
        self.weightFilePath=weightFilePath
        self.reader=ROOT.TMVA.Reader()
        self.declareDict()
        self.ready=False
        self.energyvariablelist=['lep1_pt','lep2_pt','met_et','transverse_mass_lep1_lep2','transverse_mass_lep1_met','transverse_mass_lep2_met','ptsum_lep1_lep2_met','ptsum_lep1_lep2','visible_mass','averaged_mass']

    def useRelativeVariables(self,Energy_scale):
        self.withRelativeVariables=True
        self.Energy_scale=Energy_scale
        
    def useCorrectionAsTarget(self,Corrected_variable):
	self.useCorrection=True
	self.correctionFactor=Corrected_variable

    def prepareReader(self):
        if not(self.ready):
            for varName, var in sorted(self.variables.iteritems()):
                if self.withRelativeVariables:
                    if varName==self.Energy_scale:
                        continue
                self.reader.AddVariable(varName,var[4])
            print '=================Booking Reader================='
            self.reader.BookMVA('BRT_HiggsMass',self.weightFilePath)
            print '================================================'
            self.ready=True
    def declareDict(self):
        self.variables = {}
        self.variables['lep1_pt']                      = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,200000., 40]
        self.variables['lep1_eta']                     = [''   ,'F',-10.00,10.00,  array.array('f',[0]),ROOT.TBranch(), -10.00,10.00,  200]
        self.variables['lep2_pt']                      = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,200000., 40]
        self.variables['lep2_eta']                     = [''   ,'F',-10.00,10.00,  array.array('f',[0]),ROOT.TBranch(), -10.00,10.00,  200]
        self.variables['met_et']                       = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,200000., 40]
        self.variables['transverse_mass_lep1_lep2']    = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
        self.variables['transverse_mass_lep1_met']     = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
        self.variables['transverse_mass_lep2_met']     = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
        self.variables['dphi_lep1_met']                = ['rad','F',  0.00, 3.15,  array.array('f',[0]),ROOT.TBranch(),   0.00, 3.15,   30]
        self.variables['dphi_lep2_met']                = ['rad','F',  0.00, 3.15,  array.array('f',[0]),ROOT.TBranch(),   0.00, 3.15,   30]
        self.variables['dphi_lep_lep']                 = ['rad','F',  0.00, 3.15,  array.array('f',[0]),ROOT.TBranch(),   0.00, 3.15,   30]
        self.variables['deta_lep_lep']                 = ['',   'F',  0.00,20.00,  array.array('f',[0]),ROOT.TBranch(),   0.00,20.00,  200]
        self.variables['dR_lep_lep']                   = ['',   'F',  0.00,25.00,  array.array('f',[0]),ROOT.TBranch(),   0.00,20.00,  200]
        self.variables['ptsum_lep1_lep2_met']          = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,300000., 30]
        self.variables['ptsum_lep1_lep2']              = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 25]
        self.variables['pttot_lep1_lep2_met']          = ['',   'F',  0.00, 2.00,  array.array('f',[0]),ROOT.TBranch(),   0.00, 1.10,   22]
        self.variables['pttot_lep1_lep2']              = ['',   'F',  0.00, 2.00,  array.array('f',[0]),ROOT.TBranch(),   0.00, 1.10,   22]
        self.variables['ptdiff_lep1_lep2']             = ['',   'F',  0.00, 2.00,  array.array('f',[0]),ROOT.TBranch(),   0.00, 1.10,   22]
        #self.variables['visible_mass']                 = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
        #self.variables['collinear_mass']                 = ['MeV','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   0.00,250000., 50]
        #self.variables['met_phi_centrality_smirnov_par_perp']                 = ['','F',  0.00,9999999,array.array('f',[0]),ROOT.TBranch(),   -1.45,1.45, 50]

    def mass_BRT(self,Tau1,Tau2,met_px,met_py):
        self.prepareReader()
        self.UpdateVariables(Tau1,Tau2,met_px,met_py)
	if self.useCorrection:
	    return self.reader.EvaluateMVA("BRT_HiggsMass")*variables.evalVariable(self.correctionFactor,Tau1,Tau2,met_px,met_py)/1000.
	else:
	    return self.reader.EvaluateMVA("BRT_HiggsMass")/1000.

    def UpdateVariables(self,Tau1,Tau2,met_px,met_py):
        if self.withRelativeVariables:
            self.energyscalefactor=variables.evalVariable(self.Energy_scale,Tau1,Tau2,met_px,met_py)
        for varName, var in sorted(self.variables.iteritems()):
            if varName in self.energyvariablelist and self.withRelativeVariables:
                var[4][0] = variables.evalVariable(varName,Tau1,Tau2,met_px,met_py)/self.energyscalefactor
            else:
                var[4][0] = variables.evalVariable(varName,Tau1,Tau2,met_px,met_py)
