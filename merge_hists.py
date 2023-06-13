import ROOT
import os 


path_1 = "/afs/cern.ch/user/m/msanjran/public/nominal_syst_event_bjet_og/"
path_2 = "/afs/cern.ch/user/p/pinkaew/public/nominal/"
path_out = "./systs/"
if not os.path.exists(path_out):
    os.makedirs(path_out)


files_1 = os.listdir(path_1)
files_2 = os.listdir(path_2)

for file_1 in files_1:

    for file_2 in files_2:

        if not file_1.replace(".root", "") in file_2:
            continue 
        

        # Create a new output file
        output_file = ROOT.TFile(path_out + file_1 + ".root", "RECREATE")



        # Open the input files
        file1 = ROOT.TFile.Open(path_1 + file_1, "READ")
        file2 = ROOT.TFile.Open(path_2 + file_2, "READ")

        for key in file1.GetListOfKeys():
            print "key1", key.GetName()
            if (not "bjeten" in key.GetName()) or ("bjetenls" in key.GetName()):
                continue

            hist_file1 = file1.Get(key)
            
            
            
        for key in file2.GetListOfKeys():

            print "key2", key.GetName()



            # hist_1 = file_1.Get(key)

        # Get the histograms from the input files
        histA_file1 = file1.Get("histA")
        histB_file1 = file1.Get("histB")
        histA_file2 = file2.Get("histA")
        histC_file2 = file2.Get("histC")


        # Write copies of the histograms to the output file
        output_file.cd()
        histA_file1.Write("histA")
        histB_file1.Write("histB")
        histA_file2.Write("histC")

        # Close the files
        file1.Close()
        file2.Close()
        output_file.Close()