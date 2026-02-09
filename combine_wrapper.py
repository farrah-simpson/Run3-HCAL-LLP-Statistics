import ROOT
import sys
import re
import subprocess
import argparse
import os 
import json 

ROOT.gROOT.SetBatch(True)

cwd = os.getcwd()
default_template_datacard = os.path.join( cwd, "templates/v0/datacard_TEMPLATE.txt" )

# Lifetimes (ctau) in in mm -- points to dynamically reweight to
lifetimes    = ["1000"] #["10", "30", "50", "100", "200", "300", "500", "800", "1000", "2000", "3000", "5000", "10000"]

# Temporary scale factor, otherwise get weird results
SF_temp = 0.01

# Use this for text scraping (expert only)
expected_percent = [" 2.5", "16.0", "50.0", "84.0", "97.5"] 

# ------------------------------------------------------------------------------
def parseArgs():
    """ Parse command-line arguments
    """
    parser = argparse.ArgumentParser(
        add_help=True,
        description=''
    )

    # General
    parser.add_argument("-d", "--debug",      action="store_true", default=False, help="Debug mode")
    parser.add_argument("-i", "--input",      action="store", help="Input signal file (ROOT Minituple)", required=True)
    parser.add_argument("-t", "--template",   action="store", default=default_template_datacard, help="Input template datacard")
    parser.add_argument("-f", "--filetag",    action="store", help="Input file tag", required=True)
    parser.add_argument("-c", "--ctau",       action="store", help="Input file lifetime", required=True)
    parser.add_argument("-o", "--output-dir", action="store", default="output", help="Output directory")
    parser.add_argument("--incl-score",       action="store", default=0.9, help="Signal region inclusive score cut")
    parser.add_argument("--depth-score",      action="store", default=0.8, help="Signal region depth score cut")

    args = parser.parse_args()

    return args

# ------------------------------------------------------------------------------
def calculate_bkg_prediction(tree_data_skim, lumi_sf, incl_score_cut, depth_score_cut, additional_cut_jet0="", additional_cut_jet1=""): 

    # This is really just an approximation, should use Gillian's more robust approach eventually

    if additional_cut_jet0 != "": additional_cut_jet0 = "("+additional_cut_jet0+") && "
    if additional_cut_jet1 != "": additional_cut_jet1 = "("+additional_cut_jet1+") && "

    nevents_bkg_ljdc_cr_temp       = tree_data_skim.GetEntries( additional_cut_jet0 + "(jet0_DepthTagCand == 1 && jet1_InclTagCand == 1 && jet1_scores_inc_train80 < 0.2)".format( incl_score_cut ) )
    nevents_bkg_sjdc_cr_temp       = tree_data_skim.GetEntries( additional_cut_jet1 + "(jet1_DepthTagCand == 1 && jet0_InclTagCand == 1 && jet0_scores_inc_train80 < 0.2)".format( incl_score_cut ) )
    nevents_bkg_ljdc_cr_depth_temp = tree_data_skim.GetEntries( additional_cut_jet0 + "(jet0_DepthTagCand == 1 && jet1_InclTagCand == 1 && jet0_scores_depth_LLPanywhere > {0} && jet1_scores_inc_train80 < 0.2)".format( depth_score_cut, incl_score_cut ) )
    nevents_bkg_sjdc_cr_depth_temp = tree_data_skim.GetEntries( additional_cut_jet1 + "(jet1_DepthTagCand == 1 && jet0_InclTagCand == 1 && jet1_scores_depth_LLPanywhere > {0} && jet0_scores_inc_train80 < 0.2)".format( depth_score_cut, incl_score_cut ) )
    nevents_bkg_ljdc_sr_temp       = tree_data_skim.GetEntries( additional_cut_jet0 + "(jet0_DepthTagCand == 1 && jet1_InclTagCand == 1 && jet1_scores_inc_train80 > {0} )".format( incl_score_cut ) )
    nevents_bkg_sjdc_sr_temp       = tree_data_skim.GetEntries( additional_cut_jet1 + "(jet1_DepthTagCand == 1 && jet0_InclTagCand == 1 && jet0_scores_inc_train80 > {0} )".format( incl_score_cut ) )

    """
    print( "nevents_bkg_ljdc_cr_temp", nevents_bkg_ljdc_cr_temp )       
    print( "nevents_bkg_sjdc_cr_temp", nevents_bkg_sjdc_cr_temp )        
    print( "nevents_bkg_ljdc_cr_depth_temp", nevents_bkg_ljdc_cr_depth_temp ) 
    print( "nevents_bkg_sjdc_cr_depth_temp", nevents_bkg_sjdc_cr_depth_temp ) 
    print( "nevents_bkg_ljdc_sr_temp", nevents_bkg_ljdc_sr_temp ) 
    print( "nevents_bkg_sjdc_sr_temp",nevents_bkg_sjdc_sr_temp  ) 
    """

    nevents_bkg_ljdc_srpred = lumi_sf * nevents_bkg_ljdc_sr_temp * (nevents_bkg_ljdc_cr_depth_temp / nevents_bkg_ljdc_cr_temp)
    nevents_bkg_sjdc_srpred = lumi_sf * nevents_bkg_sjdc_sr_temp * (nevents_bkg_sjdc_cr_depth_temp / nevents_bkg_sjdc_cr_temp)

    print( "nevents_bkg_ljdc_srpred",nevents_bkg_ljdc_srpred  ) 
    print( "nevents_bkg_sjdc_srpred",nevents_bkg_sjdc_srpred  ) 

    return nevents_bkg_ljdc_srpred, nevents_bkg_sjdc_srpred

def make_chain(tree_name, file_list):
    chain = ROOT.TChain(tree_name)
    for f in file_list:
        chain.Add(f)
    return chain
# ------------------------------------------------------------------------------
def main():

    # ----- Process Inputs ----- #

    args = parseArgs()

    debug = args.debug

    template_datacard = args.template
    filetag           = args.filetag
    infilepath        = args.input
    ctau_sample       = args.ctau
    incl_score_cut    = args.incl_score
    depth_score_cut   = args.depth_score
    output_dir        = args.output_dir

    unique_filetag = "{0}_{1}_{2}".format( filetag, incl_score_cut, depth_score_cut)

    ctaus        = []
    limits_obs   = []
    nevents_sig_ljdc = []
    nevents_sig_sjdc = []

    limits_expected = {}
    for val in expected_percent: limits_expected[val] = []

    # ----- Read in Data (Background Prediction) ----- #

    # TODO: Eventually should use Gillian's background estimation

    print("Reading in data tree... (this may take a few minutes)")

    # if using a partial dataset, how much to scale this up by
    lumi_sf_23 = 1.0 # 2023 all eras
    lumi_sf_22 = 6.8 # 2023D

    # currently only using a partial dataset (2023D corresponds to lumi scale factor above)
    infile_data_23 = [
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v4.1/minituple_LLPskim_2023Cv1_allscores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v4.1/minituple_LLPskim_2023Cv2_allscores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v4.1/minituple_LLPskim_2023Cv3_allscores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v4.1/minituple_LLPskim_2023Cv4_allscores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v4.1/minituple_LLPskim_2023D_allscores.root",
    ]

    infile_data_22 = ROOT.TFile.Open("/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v4.1/minituple_LLPskim_2023D_allscores.root")

    tree_data_23 = make_chain("NoSel", infile_data_23) 
    tree_data_22   = infile_data_22.Get("NoSel")

    # Copy tree but only copy these branches
    tree_data_23.SetBranchStatus("*", 0) 
    tree_data_23.SetBranchStatus("Pass_PreSel", 1) 
    tree_data_23.SetBranchStatus("jet*_DepthTagCand", 1) 
    tree_data_23.SetBranchStatus("jet*_InclTagCand", 1) 
    tree_data_23.SetBranchStatus("jet*_scores*", 1) 
    tree_data_23.SetBranchStatus("jet*_DeepCSV*", 1)

    tree_data_22.SetBranchStatus("*", 0) 
    tree_data_22.SetBranchStatus("Pass_PreSel", 1) 
    tree_data_22.SetBranchStatus("jet*_DepthTagCand", 1) 
    tree_data_22.SetBranchStatus("jet*_InclTagCand", 1) 
    tree_data_22.SetBranchStatus("jet*_scores*", 1) 
    tree_data_22.SetBranchStatus("jet*_DeepCSV*", 1)

    outfile_temp = ROOT.TFile("skim_temp_{0}.root".format(unique_filetag),"RECREATE")
    outfile_temp.cd()
    tree_data_skim_23 = tree_data_23.CopyTree("Pass_PreSel == 1")
    tree_data_skim_22 = tree_data_22.CopyTree("Pass_PreSel == 1")

    # Btag score corresponds to 2023 post BPIX (TODO: Fix)
    nevents_bkg_ljdc_srpred_23, nevents_bkg_sjdc_srpred_23               = calculate_bkg_prediction(tree_data_skim_23, lumi_sf_23, incl_score_cut, depth_score_cut) #, additional_cut_jet0="", additional_cut_jet1="")
    nevents_bkg_ljdc_srpred_btag_23, nevents_bkg_sjdc_srpred_btag_23     = calculate_bkg_prediction(tree_data_skim_23, lumi_sf_23, incl_score_cut, depth_score_cut, "jet0_DeepCSV_prob_b > 0.2435", "jet1_DeepCSV_prob_b > 0.2435") #, additional_cut_jet0="", additional_cut_jet1="")
    nevents_bkg_ljdc_srpred_nobtag_23, nevents_bkg_sjdc_srpred_nobtag_23 = calculate_bkg_prediction(tree_data_skim_23, lumi_sf_23, incl_score_cut, depth_score_cut, "jet0_DeepCSV_prob_b < 0.2435", "jet1_DeepCSV_prob_b < 0.2435") #, additional_cut_jet0="", additional_cut_jet1="")

    nevents_bkg_ljdc_srpred_22, nevents_bkg_sjdc_srpred_22               = calculate_bkg_prediction(tree_data_skim_22, lumi_sf_22, incl_score_cut, depth_score_cut) #, additional_cut_jet0="", additional_cut_jet1="")
    nevents_bkg_ljdc_srpred_btag_22, nevents_bkg_sjdc_srpred_btag_22     = calculate_bkg_prediction(tree_data_skim_22, lumi_sf_22, incl_score_cut, depth_score_cut, "jet0_DeepCSV_prob_b > 0.2435", "jet1_DeepCSV_prob_b > 0.2435") #, additional_cut_jet0="", additional_cut_jet1="")
    nevents_bkg_ljdc_srpred_nobtag_22, nevents_bkg_sjdc_srpred_nobtag_22 = calculate_bkg_prediction(tree_data_skim_22, lumi_sf_22, incl_score_cut, depth_score_cut, "jet0_DeepCSV_prob_b < 0.2435", "jet1_DeepCSV_prob_b < 0.2435") #, additional_cut_jet0="", additional_cut_jet1="")


    # ----- Read in Signal ----- #

    print("Reading in signal tree...")

    infile_sig = ROOT.TFile.Open(infilepath)
    tree_sig   = infile_sig.Get("NoSel")

    tree_sig.SetBranchStatus("*", 0) 
    tree_sig.SetBranchStatus("Pass_PreSel", 1) 
    tree_sig.SetBranchStatus("jet*_DepthTagCand", 1) 
    tree_sig.SetBranchStatus("jet*_InclTagCand", 1) 
    tree_sig.SetBranchStatus("jet*_scores*", 1) 
    tree_sig.SetBranchStatus("weight*", 1) 
    tree_sig.SetBranchStatus("LLP*", 1) 

    tree_sig_skim = tree_sig.CopyTree("Pass_PreSel == 1")

    # ----- Loop over Signal Lifetimes ----- #

    print( "Getting Event Counts...")

    for ctau_target in lifetimes:
        print( "\nCTau Target:", ctau_target )

        # Lifetime Reweight

        reweight_llp0 = "pow ( {0} / {1}, 1 ) * exp( -LLP0_DecayCtau * 10. * ( 1.0/{2} - 1.0/{3} ) )".format(ctau_sample, ctau_target, ctau_target, ctau_sample)
        reweight_llp1 = "pow ( {0} / {1}, 1 ) * exp( -LLP1_DecayCtau * 10. * ( 1.0/{2} - 1.0/{3} ) )".format(ctau_sample, ctau_target, ctau_target, ctau_sample)
        reweight = "(weight * {0} * {1})".format(reweight_llp0, reweight_llp1)

        hist_sig_ljdc = ROOT.TH1F("hist_sig_ljdc_"+ctau_target, "", 1, 0, 1)
        hist_sig_sjdc = ROOT.TH1F("hist_sig_sjdc_"+ctau_target, "", 1, 0, 1)

        tree_sig_skim.Draw("0.5 >> hist_sig_ljdc_"+ctau_target, " {0} * (jet0_DepthTagCand == 1 && jet1_InclTagCand == 1 && jet0_scores_depth_LLPanywhere > {1} && jet1_scores_inc_train80 > {2})".format(reweight, depth_score_cut, incl_score_cut ) )
        tree_sig_skim.Draw("0.5 >> hist_sig_sjdc_"+ctau_target, " {0} * (jet1_DepthTagCand == 1 && jet0_InclTagCand == 1 && jet1_scores_depth_LLPanywhere > {1} && jet0_scores_inc_train80 > {2})".format(reweight, depth_score_cut, incl_score_cut ) )

        nevents_sig_ljdc_temp = hist_sig_ljdc.Integral() * SF_temp * 100. # 100 to convert from minituple % --> net fraction 
        nevents_sig_sjdc_temp = hist_sig_sjdc.Integral() * SF_temp * 100. # 100 to convert from minituple % --> net fraction

        nevents_sig_ljdc.append( nevents_sig_ljdc_temp / SF_temp )
        nevents_sig_sjdc.append( nevents_sig_sjdc_temp / SF_temp )

        # Replace test in template datacard

        output_file = template_datacard.replace("TEMPLATE", unique_filetag + "__" + ctau_target )
        print( "Nsig (ljdc, sjdc):", nevents_sig_ljdc_temp, nevents_sig_sjdc_temp)

        replacements = {
            "SIGLJDC": f"{nevents_sig_ljdc_temp:04.2f}", 
            "SIGSJDC": f"{nevents_sig_sjdc_temp:04.2f}",
            "BKGLJDC_23_B": f"{nevents_bkg_ljdc_srpred_btag_23:04.2f}",
            "BKGLJDC_23_XB": f"{nevents_bkg_ljdc_srpred_nobtag_23:04.2f}",
            "BKGSJDC_23_B": f"{nevents_bkg_sjdc_srpred_btag_23:04.2f}",
            "BKGSJDC_23_XB": f"{nevents_bkg_sjdc_srpred_nobtag_23:04.2f}",
            "BKGLJDC_23": f"{nevents_bkg_ljdc_srpred_23:04.2f}", 
            "BKGSJDC_23": f"{nevents_bkg_sjdc_srpred_23:04.2f}",
            "BKGLJDC_22_B": f"{nevents_bkg_ljdc_srpred_btag_22:04.2f}",
            "BKGLJDC_22_XB": f"{nevents_bkg_ljdc_srpred_nobtag_22:04.2f}",
            "BKGSJDC_22_B": f"{nevents_bkg_sjdc_srpred_btag_22:04.2f}",
            "BKGSJDC_22_XB": f"{nevents_bkg_sjdc_srpred_nobtag_22:04.2f}"
            "BKGLJDC_22": f"{nevents_bkg_ljdc_srpred_22:04.2f}", 
            "BKGSJDC_22": f"{nevents_bkg_sjdc_srpred_22:04.2f}",

        }

        pattern = re.compile("|".join(re.escape(k) for k in replacements))

        with open(template_datacard) as fin, open(output_file, "w") as fout:
            for line in fin:
                fout.write(pattern.sub(lambda m: replacements[m.group(0)], line))

        output = subprocess.check_output("combine -M AsymptoticLimits {}".format(output_file), shell=True, text=True)
        #print( output )

        match_all = True 
       
        match = re.search(r"Observed Limit:\s*r\s*<\s*([0-9.]+)", output)

        if match:
            ctaus.append( float(ctau_target) )
            limits_obs.append( float(match.group(1)) * SF_temp )
        else:
            limits_obs.append( -1 )
            match_all = False

        for val in expected_percent: 
            pattern = rf"Expected {val}%:\s*r\s*<\s*([0-9.]+)"
            match = re.search(pattern, output)
            if match:
                limits_expected[val].append( float( match.group(1) ) * SF_temp )
            else:
                limits_expected[val].append( -1 )
                match_all = False

        if not match_all: 
            print("WARNING: could not extract all limit information for:", ctau_target, "(more info available in debug mode)" )
            if debug: print( output )

    data = {}
    data["ctaus"] = ctaus
    data["limits_obs"] = limits_obs
    data["limits_exp"] = limits_expected
    data["nevents_sig_ljdc"] = nevents_sig_ljdc
    data["nevents_sig_sjdc"] = nevents_sig_sjdc
    data["nevents_bkg_ljdc_23"] = nevents_bkg_ljdc_srpred_23
    data["nevents_bkg_sjdc_23"] = nevents_bkg_sjdc_srpred_23
    data["nevents_bkg_ljdc_22"] = nevents_bkg_ljdc_srpred_22
    data["nevents_bkg_sjdc_22"] = nevents_bkg_sjdc_srpred_22


    if not os.path.exists(output_dir): 
        os.makedirs(output_dir)

    outfile_path = os.path.join( output_dir, "{0}_inc{1}_depth{2}.json".format(filetag, incl_score_cut, depth_score_cut ) )

    with open(outfile_path, "w") as f:
        json.dump(data, f, indent=2)

    print( "--------------------------------------" )
    print( "CTaus: ", ctaus )
    print( "Limits:", limits_expected["50.0"] )
    print( "LJDC:  ", nevents_sig_ljdc )
    print( "SJDC:  ", nevents_sig_sjdc )
    print( "--------------------------------------" )
    print( "Json file written to:", outfile_path )

if __name__ == '__main__':
    main()
