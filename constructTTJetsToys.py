import ROOT 
import numpy as np


path_MC_TTJets= "/afs/cern.ch/user/p/phnattla/public/Uncertainties/TTJets.root"
path_MC_TTJets = "./RebinnedHists_fixed/TTJets.root"
file_MC = ROOT.TFile.Open(path_MC_TTJets, "READ")


path_MC_bkg= "/afs/cern.ch/user/p/phnattla/public/Uncertainties/bkg_merged.root"
path_MC_bkg = "./RebinnedHists_fixed/bkg.root"
file_MC_bkg = ROOT.TFile.Open(path_MC_bkg, "READ")
hist_MC_bkg = file_MC_bkg.Get("bjeten_nominal")

hist_nom = file_MC.Get("bjeten_nominal")
hist_up = file_MC.Get("bjeten_mt_Up")
hist_down = file_MC.Get("bjeten_mt_Down")

n_bins = hist_nom.GetNbinsX()
x_min = hist_nom.GetXaxis().GetXmin()
x_max = hist_nom.GetXaxis().GetXmax()

toy_hist_nom = ROOT.TH1F("bjeten_nominal", "Toy Histogram (Nominal)", n_bins, x_min, x_max)
toy_hist_up = ROOT.TH1F("bjeten_mt_Up", "Toy Histogram (Up)", n_bins, x_min, x_max)
toy_hist_down = ROOT.TH1F("bjeten_mt_Down", "Toy Histogram (Down)", n_bins, x_min, x_max)






for i in range(hist_nom.GetNbinsX()):

    bkg = hist_MC_bkg.GetBinContent(i+1)

    nom = hist_nom.GetBinContent(i+1)
    up = hist_up.GetBinContent(i+1)
    down = hist_down.GetBinContent(i+1)

    toy_nom = np.random.poisson(nom+bkg)
    toy_up = np.random.poisson(up+bkg)
    toy_down = np.random.poisson(down+bkg)

    toy_hist_nom.SetBinContent(i+1, toy_nom)
    toy_hist_up.SetBinContent(i+1, toy_up)
    toy_hist_down.SetBinContent(i+1, toy_down)


output_file = ROOT.TFile("data_obs.root", "RECREATE")
output_file.cd()
toy_hist_nom.Write()
toy_hist_up.Write()
toy_hist_down.Write()
output_file.Close()
