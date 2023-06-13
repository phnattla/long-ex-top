import ROOT, os, array


path = "/afs/cern.ch/user/p/pinkaew/public/systematics/"
outpath = "./RebinnedHists/"
outpath_fixed = "./RebinnedHists_fixed/"
if not os.path.exists(outpath):
    os.makedirs(outpath)
if not os.path.exists(outpath_fixed):
    os.makedirs(outpath_fixed)
processes = ["TTJets", "WW", "tW", "atW", "DY", "MuonEG_2016"]


def fixUncertainties(h, hUp, hDown):
    hUpOut = hUp.Clone()
    hDownOut = hDown.Clone()

    for ibin in range(1, hUpOut.GetNbinsX()+1):
        upvarContent = hUp.GetBinContent(ibin)
        downvarContent = hDown.GetBinContent(ibin)
        er = (upvarContent-downvarContent)/2.
        nomVal = h.GetBinContent(ibin)
        hUpOut.SetBinContent(ibin, nomVal+er)
        hDownOut.SetBinContent(ibin, nomVal-er)
        hUpOut.SetDirectory(0)
        hDownOut.SetDirectory(0)
    return hUpOut,hDownOut





def frangeN(start, end, N):
    a = []
    tmp = start
    step = abs(end - start) / N
    while(tmp <= end):
        a.append(tmp)
        tmp += step
    return a


def handle_overflow_underflow(hist):

    n_bins = hist.GetNbinsX()

    for bin in range(1, n_bins+1):
        if hist.GetBinContent(bin) < 1e-3:
            hist.SetBinContent(bin, 1e-3)
            hist.SetBinError(bin, 1e-3)

    underflow = hist.GetBinContent(0)
    overflow = hist.GetBinContent(n_bins+1)

    hist.SetBinContent(1, hist.GetBinContent(1)+underflow)
    hist.SetBinContent(n_bins, hist.GetBinContent(n_bins)+overflow)

    hist.SetBinContent(0, 0)
    hist.SetBinContent(n_bins+1, 0)
    # bin error 

    return hist 


for process in processes:

    # Open the input file
    input_file = ROOT.TFile.Open(path + process + ".root", "READ")
    output_file = ROOT.TFile(outpath+process+".root", "RECREATE")

    for hist_key in input_file.GetListOfKeys():

        if "bjetenls" in hist_key.GetName(): continue

        print hist_key.GetName()

    
        # Get the input histogram
        input_hist = input_file.Get(hist_key.GetName())

        # Define the rebinning parameters
        start_bin = 4    # Exclude the first two bins

        # Get the number of bins and adjust for rebinning
        n_bins = input_hist.GetNbinsX()
        rebin_n_bins = int((n_bins - start_bin + 1) )

        # Create the output histogram with the adjusted binning
        # output_hist = ROOT.TH1F(hist_key.GetName(), "Rebinned Histogram", rebin_n_bins, 30, 300)
        output_hist = input_hist.Clone()

        # Loop over the bins of the input histogram, excluding the first two bins
        # for i in range(start_bin, n_bins + 1):
        #     bin_content = input_hist.GetBinContent(i)
        #     bin_error = input_hist.GetBinError(i)
        #     rebin_index = int((i - start_bin) ) + 1
        #     output_hist.SetBinContent(rebin_index, bin_content)
        #     output_hist.SetBinError(rebin_index, bin_error)

        bins = frangeN(30, 300, 9) #[30*i for i in range(1, n_bins+1)]
        #bins = [30,60,90,120,150,300]
        bins = array.array("d", bins)

        # print(bins); exit()
        output_hist = output_hist.Rebin(len(bins)-1, hist_key.GetName(), bins)

        output_hist = handle_overflow_underflow(output_hist)

        # Create the output file and write the rebinned histogram
        
        output_file.cd()
        output_hist.Write()
input_file.Close()
for process in processes:

    # Open the input file
    input_file = ROOT.TFile.Open(outpath + process + ".root", "READ")
    output_file = ROOT.TFile(outpath_fixed+process+".root", "RECREATE")
    for histkey in input_file.GetListOfKeys():
        if "bjetenls" in histkey.GetName(): 
            continue
        if ("Up" in histkey.GetName()) and ("fsr" in histkey.GetName() or "isr" in histkey.GetName() or "scale" in histkey.GetName() or "mt" in histkey.GetName() or "jec" in histkey.GetName()):
            name = histkey.GetName()
            inputUp = input_file.Get(name)
            inputDown = input_file.Get(name.replace("_Up","_Down"))
            inputUp.SetDirectory(0)
            inputDown.SetDirectory(0)

            inputNom = input_file.Get("bjeten_nominal")
            upout, downout = fixUncertainties(inputNom, inputUp, inputDown)
            output_file.cd()
            upout.Write()
            downout.Write()
        else:
            name = histkey.GetName()
            input = input_file.Get(name)
            input.SetDirectory(0)
            output_file.cd()
            input.Write()

        # Close the files
    output_file.Close()

input_file.Close()
