import ROOT
import sys
import re
import subprocess
import argparse
import os 
import json 
import tempfile
import glob

ROOT.gROOT.SetBatch(True)

cwd = os.getcwd()
default_template_datacard = os.path.join( cwd, "templates/v1/datacard_TEMPLATE.txt" )

# Lifetimes (ctau) in in mm -- points to dynamically reweight to
lifetimes    = ["10", "30", "50", "100", "200", "300", "500", "800", "1000", "2000", "3000", "5000", "10000"]

# Temporary scale factor applied to signal yield to stabilize results in combine
#limits are rescale back, so final results are unaffected otherwise get weird results
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
    parser.add_argument("--input-jer-up", action="store", help="Input signal ROOT file with JER up", required=True)
    parser.add_argument("--input-jer-down", action="store", help="Input signal ROOT file with JER down", required=True)

    args = parser.parse_args()

    return args

# ------------------------------------------------------------------------------
### INSERT THE VALUES FROM TABLE BELOW FOR NOW
def calculate_bkg_prediction(tree_data_skim, lumi_sf, incl_score_cut, depth_score_cut, additional_cut_jet0="", additional_cut_jet1=""):

    if additional_cut_jet0 != "": additional_cut_jet0 = f"({additional_cut_jet0}) && "
    if additional_cut_jet1 != "": additional_cut_jet1 = f"({additional_cut_jet1}) && "

    nevents_bkg_ljdc_cr_temp       = tree_data_skim.GetEntries(additional_cut_jet0 + "(jet0_DepthTagCand == 1 && jet1_InclTagCand == 1 && jet1_scores_inc_train80 < 0.2)")
    nevents_bkg_sjdc_cr_temp       = tree_data_skim.GetEntries(additional_cut_jet1 + "(jet1_DepthTagCand == 1 && jet0_InclTagCand == 1 && jet0_scores_inc_train80 < 0.2)")
    nevents_bkg_ljdc_cr_depth_temp = tree_data_skim.GetEntries(additional_cut_jet0 + "(jet0_DepthTagCand == 1 && jet1_InclTagCand == 1 && jet0_scores_depth_LLPanywhere > {0} && jet1_scores_inc_train80 < 0.2)".format(depth_score_cut))
    nevents_bkg_sjdc_cr_depth_temp = tree_data_skim.GetEntries(additional_cut_jet1 + "(jet1_DepthTagCand == 1 && jet0_InclTagCand == 1 && jet1_scores_depth_LLPanywhere > {0} && jet0_scores_inc_train80 < 0.2)".format(depth_score_cut))
    nevents_bkg_ljdc_sr_temp       = tree_data_skim.GetEntries(additional_cut_jet0 + "(jet0_DepthTagCand == 1 && jet1_InclTagCand == 1 && jet1_scores_inc_train80 > {0})".format(incl_score_cut))
    nevents_bkg_sjdc_sr_temp       = tree_data_skim.GetEntries(additional_cut_jet1 + "(jet1_DepthTagCand == 1 && jet0_InclTagCand == 1 && jet0_scores_inc_train80 > {0})".format(incl_score_cut))

    lj_ratio = (nevents_bkg_ljdc_cr_depth_temp / nevents_bkg_ljdc_cr_temp) if nevents_bkg_ljdc_cr_temp else 0.0
    sj_ratio = (nevents_bkg_sjdc_cr_depth_temp / nevents_bkg_sjdc_cr_temp) if nevents_bkg_sjdc_cr_temp else 0.0

    nevents_bkg_ljdc_srpred = lumi_sf * nevents_bkg_ljdc_sr_temp * lj_ratio
    nevents_bkg_sjdc_srpred = lumi_sf * nevents_bkg_sjdc_sr_temp * sj_ratio

    return nevents_bkg_ljdc_srpred, nevents_bkg_sjdc_srpred

def make_chain(tree_name, file_list):
    chain = ROOT.TChain(tree_name)
    for f in file_list:
        chain.Add(f)
    return chain

def get_signal_yield(infilepath, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc, pileupweight_="nominal", sys="nominal"):

    # if using a partial dataset, how much to scale this up by
    lumi_sf_23 = 1.0 
    lumi_sf_22 = 1.0 
    
    lumi_2022 = 31.51
    lumi_2023 =28.91
    lumi_total = lumi_2022 + lumi_2023
 
    # ----- Read in Signal ----- #
 
    print("Reading in signal tree...")
 
    infile_sig = ROOT.TFile.Open(infilepath) 
    tree_sig  = infile_sig.Get("NoSel")
 
    tree_sig.SetBranchStatus("*", 0) 
    tree_sig.SetBranchStatus("Pass_PreSel", 1) 
    tree_sig.SetBranchStatus("jet*_DepthTagCand", 1) 
    tree_sig.SetBranchStatus("jet*_InclTagCand", 1) 
    tree_sig.SetBranchStatus("jet*_scores*", 1) 
    tree_sig.SetBranchStatus("jet0_Pt", 1) 
    tree_sig.SetBranchStatus("jet1_Pt", 1) 
    tree_sig.SetBranchStatus("weight*", 1) 
    tree_sig.SetBranchStatus("event_weight", 1) 
    tree_sig.SetBranchStatus("LLP*", 1) 
    tree_sig.SetBranchStatus("L1_prescale_weight", 1) 
    tree_sig.SetBranchStatus("puWeight*", 1) 
    tree_sig.SetBranchStatus("jet0_jet1_dPhi", 1)
    tree_sig.SetBranchStatus("Flag_METFilters_2022_2023_PromptReco", 1)
    tree_sig.SetBranchStatus("Pass_HLTDisplacedJet", 1) 
    
    train_frac = 0.2
    lj_train_cut = "(int(jet0_Pt * 1000) % 10) >= 8"
    sj_train_cut = "(int(jet1_Pt * 1000) % 10) >= 8"
 
    deltaPhi_cut = "(abs(jet0_jet1_dPhi) > 0.2)"

    # Lifetime Reweight

    reweight_llp0 = "pow ( {0} / {1}, 1 ) * exp( -LLP0_DecayCtau * 10. * ( 1.0/{2} - 1.0/{3} ) )".format(ctau_sample, ctau_target, ctau_target, ctau_sample)
    reweight_llp1 = "pow ( {0} / {1}, 1 ) * exp( -LLP1_DecayCtau * 10. * ( 1.0/{2} - 1.0/{3} ) )".format(ctau_sample, ctau_target, ctau_target, ctau_sample)
    #reweight = "( L1_prescale_weight * event_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)
    if pileupweight_ == "nominal": reweight = "( puWeight * L1_prescale_weight * event_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)
    elif pileupweight_ == "pileupWeightUp": reweight = "( puWeightUp * L1_prescale_weight * event_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)
    elif pileupweight_ == "pileupWeightDown": reweight = "( puWeightDown * L1_prescale_weight * event_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)

    hist_sig_ljdc = ROOT.TH1F(
        f"hist_sig_ljdc_{ctau_target}_{sys}_{pileupweight_}",
        "",
        1,0,1
    )
    
    hist_sig_sjdc = ROOT.TH1F(
        f"hist_sig_sjdc_{ctau_target}_{sys}_{pileupweight_}",
        "",
        1,0,1
    )

    print("tree entries:", tree_sig.GetEntries())
    print("reweight =", reweight)

    tree_sig.Draw(
            "0.5 >> hist_sig_ljdc_"+ctau_target+"_"+sys+"_"+pileupweight_,
       " ({0}) * ( Pass_HLTDisplacedJet == 1 && Pass_PreSel == 1 && {1} && jet0_DepthTagCand == 1 && jet1_InclTagCand == 1 && jet0_scores_depth_LLPanywhere > {2} && jet1_scores_inc_train80 > {3} && {4})".format(
           reweight, lj_train_cut, lj_depth, lj_inc, deltaPhi_cut
       )
    )

    tree_sig.Draw(
       "0.5 >> hist_sig_sjdc_"+ctau_target+"_"+sys+"_"+pileupweight_,
       " ({0}) * ( Pass_HLTDisplacedJet == 1 && Pass_PreSel == 1 && {1} && jet1_DepthTagCand == 1 && jet0_InclTagCand == 1 && jet1_scores_depth_LLPanywhere > {2} && jet0_scores_inc_train80 > {3} && {4})".format(
           reweight, sj_train_cut, sj_depth, sj_inc, deltaPhi_cut
       )
    )

    #combined 2022+2023
    nevents_sig_ljdc_temp = hist_sig_ljdc.Integral() * SF_temp * 100. / train_frac # 100 to convert from minituple % --> net fraction 
    nevents_sig_sjdc_temp = hist_sig_sjdc.Integral() * SF_temp * 100. /train_frac # 100 to convert from minituple % --> net fraction

    infile_sig.Close()

    return {
    "SIGLJDC_22": nevents_sig_ljdc_temp * lumi_2022 / lumi_total,
    "SIGLJDC_23": nevents_sig_ljdc_temp * lumi_2023 / lumi_total,
    "SIGSJDC_22": nevents_sig_sjdc_temp * lumi_2022 / lumi_total,
    "SIGSJDC_23": nevents_sig_sjdc_temp * lumi_2023 / lumi_total,
    }


# ------------------------------------------------------------------------------
def main():

    # ----- Process Inputs ----- #
    args = parseArgs()

    debug = args.debug

    template_datacard = args.template
    filetag           = args.filetag
    infilepath        = args.input
    ctau_sample       = args.ctau
    output_dir        = args.output_dir

    bkg_table = [
        {
            "lj_depth": 0.99,
            "lj_inc": 0.98,
            "sj_depth": 0.98,
            "sj_inc": 0.9,
            "bkg": {
                2022: {"lj": 0.02, "sj": 0.09},
                2023: {"lj": 0.11, "sj": 0.38},
            },
            "bkg_err": {
                2022: {"lj": 0.01, "sj": 0.03},
                2023: {"lj": 0.02, "sj": 0.22},
            },
        },
#            "lj_depth": 0.965,
#            "lj_inc": 0.695,
#            "sj_depth": 0.965,
#            "sj_inc": 0.315,
#            "CR": 0.20,
#            "bkg": {
#                2022: {"lj": 13.56, "sj": 21.45, "comb": 37.15},
#                2023: {"lj": 0.58,  "sj": 16.56, "comb": 97.65},
#            },
#            "bkg_err": {
#                2022: {"lj": 1.40, "sj": 4.13, "comb": 3.27},
#                2023: {"lj": 0.03, "sj": 2.86, "comb": 3.50},
#            },
#        },
#        {
#            "lj_depth": 0.965,
#            "lj_inc": 0.845,
#            "sj_depth": 0.975,
#            "sj_inc": 0.375,
#            "CR": 0.20,
#            "bkg": {
#                2022: {"lj": 5.27, "sj": 16.01, "comb": 25.23},
#                2023: {"lj": 9.79, "sj": 16.56, "comb": 38.38},
#            },
#            "bkg_err": {
#                2022: {"lj": 0.54, "sj": 3.38, "comb": 2.25},
#                2023: {"lj": 0.72, "sj": 2.86, "comb": 2.51},
#            },
#        },
#        {
#            "lj_depth": 0.975,
#            "lj_inc": 0.415,
#            "sj_depth": 0.975,
#            "sj_inc": 0.415,
#            "CR": 0.20,
#            "bkg": {
#                2022: {"lj": 24.01, "sj": 14.59, "comb": 38.52},
#                2023: {"lj": 43.68, "sj": 15.04, "comb": 61.36},
#            },
#            "bkg_err": {
#                2022: {"lj": 2.90, "sj": 3.08, "comb": 3.98},
#                2023: {"lj": 3.77, "sj": 2.60, "comb": 4.72},
#            },
#        },
#        {
#            "lj_depth": 0.95,
#            "lj_inc": 0.97,
#            "sj_depth": 0.95,
#            "sj_inc": 0.97,
#            "CR": 0.20,
#            "bkg": {
#                2022: {"lj": 1.00, "sj": 0.36, "comb": 1.38},
#                2023: {"lj": 1.65,  "sj": 0.71, "comb": 2.33},
#            },
#            "bkg_err": {
#                2022: {"lj": 0.08, "sj": 0.05, "comb": 0.1},
#                2023: {"lj": 0.11, "sj": 0.11, "comb": 0.13},
#            },
#        },
 
    ]
    print("Reading in data tree... (this may take a few minutes)")
   
    infile_data_23 = [
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Cv1_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Cv2_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Cv3_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Cv4_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Dv1_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Dv2_scores.root",
    ]
    infile_data_22 = [
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2022Dv1_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2022Ev1_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2022Fv1_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2022Gv1_scores.root",
    ]
    
    tree_data_23 = make_chain("NoSel", infile_data_23) 
    tree_data_22   = make_chain("NoSel", infile_data_22)
    
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

    for row in bkg_table:

        lj_depth = row["lj_depth"]
        lj_inc   = row["lj_inc"]
        sj_depth = row["sj_depth"]
        sj_inc   = row["sj_inc"]

        unique_filetag = "{0}_inc_{1}_{2}_depth_{3}_{4}".format( filetag, lj_inc, sj_inc, lj_depth, sj_depth)
    
        ctaus        = []
        limits_obs   = []
        nevents_sig_ljdc_23 = []
        nevents_sig_sjdc_23 = []
        nevents_sig_ljdc_22 = []
        nevents_sig_sjdc_22 = []
    
        limits_expected = {}
        for val in expected_percent: limits_expected[val] = []
 
        # ----- Read in Data (Background Prediction) ----- #

        # Keep temp ROOT file out of your quota-heavy working dir, and clean it up after use
        tmpdir = tempfile.mkdtemp(prefix="llp_limits_")
        tmp_root = os.path.join(tmpdir, f"skim_temp_{unique_filetag}.root")
        outfile_temp = ROOT.TFile(tmp_root, "RECREATE")

        outfile_temp.cd()
        tree_data_skim_23 = tree_data_23.CopyTree("Pass_PreSel == 1")
        tree_data_skim_22 = tree_data_22.CopyTree("Pass_PreSel == 1")
      
        nevents_bkg_ljdc_srpred_23, nevents_bkg_sjdc_srpred_23  = row["bkg"][2023]["lj"], row["bkg"][2023]["sj"]
        nevents_bkg_ljdc_srpred_22, nevents_bkg_sjdc_srpred_22  = row["bkg"][2022]["lj"], row["bkg"][2022]["sj"] 
     
        # ----- Loop over Signal Lifetimes ----- #
        print( "Getting Event Counts...")

        for ctau_target in lifetimes:

            print( "\nCTau Target:", ctau_target )

            nominal = get_signal_yield(args.input, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc)

            nevents_sig_ljdc_23.append( nominal["SIGLJDC_23"] / SF_temp )
            nevents_sig_sjdc_23.append( nominal["SIGSJDC_23"] / SF_temp )
            nevents_sig_ljdc_22.append( nominal["SIGLJDC_22"] / SF_temp )
            nevents_sig_sjdc_22.append( nominal["SIGSJDC_22"] / SF_temp )

            jer_up = get_signal_yield(args.input_jer_up, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc, "nominal", "jer_up")
            jer_down = get_signal_yield(args.input_jer_down, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc, "nominal", "jer_down")

            jer = {}
            for key in nominal:
                nom = nominal[key]
                up = jer_up[key]
                down = jer_down[key]

                if nom <= 0:
                    val = 1.0
                elif nom < 1e-3: #guard against very small yields
                    if debug: print( "guarding against small yield:", nom )
                    val = 1.0
                else:
                    val = max(up/nom, nom/down if down > 0 else 1.0)
                if val > 2.0: print( "WARNING LARGE JER:", val )
                MAX_jer = 2.0  # 100% uncertainty cap
                val = min(val, MAX_jer) # For ctau = 100,200 HToSSTo4B_125_50
                
                jer[key] = val

                print("nom, up, down:", nominal[key], jer_up[key], jer_down[key])

            pu_up = get_signal_yield(args.input, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc, "pileupWeightUp", "nominal")
            pu_down = get_signal_yield(args.input, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc, "pileupWeightDown", "nominal")

            pu = {}
            for pukey in nominal:
                punom = nominal[pukey]
                up = pu_up[pukey]
                down = pu_down[pukey]

                if punom <= 0:
                    val = 1.0
                elif punom < 1e-3: #guard against very small yields
                    if debug: print( "guarding against small yield:", punom )
                    val = 1.0
                else:
                    val = max(up/punom, punom/down if down > 0 else 1.0)
                if val > 2.0: print( "WARNING LARGE PU:", val )
                MAX_PU = 2.0  # 100% uncertainty cap
                val = min(val, MAX_PU) # For ctau = 100,200 HToSSTo4B_125_50
                
                pu[pukey] = val

                print("nom, pileup up, pileup down:", nominal[pukey], pu_up[pukey], pu_down[pukey])

            # Replace test in template datacard
            output_file = template_datacard.replace("TEMPLATE", unique_filetag + "__" + ctau_target )

            print("Nevents LJDC 23:", nominal["SIGLJDC_23"])
            print("Nevents SJDC 23:", nominal["SIGSJDC_23"])
            print("Nevents LJDC 22:", nominal["SIGLJDC_22"])
            print("Nevents SJDC 22:", nominal["SIGSJDC_22"])

            replacements = {
                "SIGLJDC_23": f"{nominal['SIGLJDC_23']:.6e}", 
                "SIGSJDC_23": f"{nominal['SIGSJDC_23']:.6e}",
                "SIGLJDC_22": f"{nominal['SIGLJDC_22']:.6e}", 
                "SIGSJDC_22": f"{nominal['SIGSJDC_22']:.6e}",

                "JERLJDC_23": f"{jer['SIGLJDC_23']:.6e}",
                "JERSJDC_23": f"{jer['SIGSJDC_23']:.6e}",
                "JERLJDC_22": f"{jer['SIGLJDC_22']:.6e}",
                "JERSJDC_22": f"{jer['SIGSJDC_22']:.6e}",

                "PULJDC_23": f"{pu['SIGLJDC_23']:.6e}",
                "PUSJDC_23": f"{pu['SIGSJDC_23']:.6e}",
                "PULJDC_22": f"{pu['SIGLJDC_22']:.6e}",
                "PUSJDC_22": f"{pu['SIGSJDC_22']:.6e}",

                "BKGLJDC_23": f"{nevents_bkg_ljdc_srpred_23:.6e}", 
                "BKGSJDC_23": f"{nevents_bkg_sjdc_srpred_23:.6e}",
                "BKGLJDC_22": f"{nevents_bkg_ljdc_srpred_22:.6e}", 
                "BKGSJDC_22": f"{nevents_bkg_sjdc_srpred_22:.6e}",
            }

            pattern = re.compile("|".join(re.escape(k) for k in replacements))

            with open(template_datacard) as fin, open(output_file, "w") as fout:
                for line in fin:
                    fout.write(pattern.sub(lambda m: replacements[m.group(0)], line))

            output = subprocess.check_output("combine -M AsymptoticLimits {}".format(output_file), shell=True, text=True)

            # remove combine ROOT byproducts
            for f in glob.glob("higgsCombine*.root"):
                try:
                    os.remove(f)
                except OSError:
                    pass

            match_all = True 
       
            match = re.search(r"Observed Limit:\s*r\s*<\s*([0-9.]+)", output)
            ctaus.append( float(ctau_target) )

            if match:
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

            outfile_temp.Close()
            if os.path.exists(tmp_root):
                os.remove(tmp_root)
            try:
                os.rmdir(tmpdir)
            except OSError:
                pass


        data = {}
        data["ctaus"] = ctaus
        data["limits_obs"] = limits_obs
        data["limits_exp"] = limits_expected
        data["nevents_sig_ljdc_23"] = nevents_sig_ljdc_23
        data["nevents_sig_sjdc_23"] = nevents_sig_sjdc_23
        data["nevents_sig_ljdc_22"] = nevents_sig_ljdc_22
        data["nevents_sig_sjdc_22"] = nevents_sig_sjdc_22
        data["nevents_bkg_ljdc_23"] = nevents_bkg_ljdc_srpred_23
        data["nevents_bkg_sjdc_23"] = nevents_bkg_sjdc_srpred_23
        data["nevents_bkg_ljdc_22"] = nevents_bkg_ljdc_srpred_22
        data["nevents_bkg_sjdc_22"] = nevents_bkg_sjdc_srpred_22
    
        if not os.path.exists(output_dir): 
            os.makedirs(output_dir)
    
        outfile_path = os.path.join( output_dir, "{0}_inc{1}_{2}_depth{3}_{4}.json".format(filetag, lj_inc, sj_inc, lj_depth, sj_depth ) )
    
        with open(outfile_path, "w") as f:
            json.dump(data, f, indent=2)
    
        print( "--------------------------------------" )
        print( "CTaus: ", ctaus )
        print( "Limits:", limits_expected["50.0"] )
        print( "LJDC 23:  ", nevents_sig_ljdc_23 )
        print( "SJDC 23:  ", nevents_sig_sjdc_23 )
        print( "LJDC 22:  ", nevents_sig_ljdc_22 )
        print( "SJDC 22:  ", nevents_sig_sjdc_22 )
    
        print( "--------------------------------------" )
        print( "Json file written to:", outfile_path )

if __name__ == '__main__':
    main()
