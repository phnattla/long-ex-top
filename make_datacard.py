import pdb 

def generate_datacard():
    # bins = ["b1", "b2", "b3+"]
    bins = ["n_bj_1o2"]
    processes = ["TTJets", "WW", "tW", "atW", "DY"]
    nuisances = ["lumi", "pileup", "lep_eff"] + ["jec_{0}".format(i) for i in range(27)]
    top_nuisances = ["mt", "fsr", "isr", "rf_scales", "tpt_reweight"]
    xSec_nuisances = [0, 1.5, 1.25, 1.25, 1.5]

    datacard = ""

    # Datacard header
    datacard += "imax {0}\n".format(len(bins))
    datacard += "jmax {0}\n".format(len(processes) - 1)
    datacard += "kmax *\n"
    datacard += "-" * 50 + "\n"
    for process in processes:
        # datacard += "shapes {} * /afs/cern.ch/user/p/phnattla/public/Uncertainties/{}.root bjeten_nominal bjeten_$SYSTEMATIC\n".format(process, process)
        datacard += "shapes {} * /eos/user/f/fmausolf/CMSDAS2023_TopMass/CMSSW_10_2_13/src/RebinnedHists_fixed/{}.root bjeten_nominal bjeten_$SYSTEMATIC\n".format(process, process)
    #datacard += "\nshapes data_obs * data_obs.root bjeten_nominal\n"
    datacard += "\nshapes data_obs * /eos/user/f/fmausolf/CMSDAS2023_TopMass/CMSSW_10_2_13/src/RebinnedHists/MuonEG_2016.root bjeten_nominal\n"
    datacard += "-" * 50 + "\n"

    # Observations
    datacard += "bin\t\t\t\t{0}\n".format("\t".join(bins))
    datacard += "observation\t\t\t{0}\n".format("\t".join(["-1"] * len(bins)))  # Replace -1 with actual observation values
    datacard += "-" * 50 + "\n"

    # Bin and process definitions
    datacard += "bin\t\t\t\t{0}\n".format("\t".join([bin_name for bin_name in bins for _ in processes]))
    datacard += "process\t\t\t\t{0}\n".format("\t".join(processes * len(bins)))
    datacard += "process\t\t\t\t{0}\n".format("\t".join([str(i) for i in range(len(processes))] * len(bins)))
    datacard += "rate\t\t\t\t{0}\n".format("\t".join(["-1"] * (len(processes) * len(bins))))  # Replace -1 with actual rates
    datacard += "-" * 50 + "\n"

    # Nuisance parameter definitions
    for nuisance in top_nuisances:
        if nuisance == "mt":
            datacard += "{0}_\tshapeU\t\t\t{1}".format(nuisance, "0.3333333\t-\t-\t-\t-\t"*len(bins)) + "\n"
        else:
            datacard += "{0}_\tshape\t\t\t{1}".format(nuisance, "1\t-\t-\t-\t-\t"*len(bins)) + "\n"

    for nuisance in nuisances:
        if nuisance == "lumi":
            datacard += "{0}_\tlnN\t\t\t{1}".format(nuisance, "1.024\t"*len(processes)*len(bins)) + "\n"
        else:
            datacard += "{0}_\tshape\t\t\t{1}".format(nuisance, "1\t"*len(processes)*len(bins)) + "\n"
    for proc_index, nuisance in enumerate(xSec_nuisances):
        if nuisance >0:
            datacard += "{0}_\tlnN\t\t\t{1}".format("xSec_" + processes[proc_index],
                "-\t"*proc_index +
                str(nuisance) + "\t" +
                "-\t" * (len(processes)-proc_index-1)
            ) + "\n"




    datacard += "-" * 50 + "\n"

    datacard += "* autoMCStats 10 0 1" + "\n"

    return datacard


datacard_content = generate_datacard()
file_name = "datacard_new.txt"  # Name of the output file

with open(file_name, "w") as file:
    file.write(datacard_content)

print("Datacard saved to {0}".format(file_name))

