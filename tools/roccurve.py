import ROOT
#Make roc curve over limited range
def roccurve(hist1,hist2):
    nBins=hist1.GetNbinsX()
    if not(nBins==hist2.GetNbinsX()):
        print "Error, uneven histogram bin numbers."
        return
    hist1_temp=hist1.Clone('hist1')
    hist2_temp=hist2.Clone('hist2')
    hist1_temp.Scale(1./hist1_temp.Integral())
    hist2_temp.Scale(1./hist2_temp.Integral())
    Roc=ROOT.TGraph(nBins+2)
    hist1_passing=1.
    hist2_passing=1.

    for i in xrange(nBins+2):
        Roc.SetPoint(i,hist1_passing,1.-hist2_passing)
        hist1_passing-=hist1_temp.GetBinContent(i)
        hist2_passing-=hist2_temp.GetBinContent(i)
    return Roc
