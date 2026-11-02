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

def get_signal_yield(infilepath, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc, pileupweight_="nominal", sys="nominal",trig_cut="1", require_hlt_decision=True):

    # scale the signal based on the difference in signal efficiencies in pre and post BPix
    sigeff_sf_23 = 0.72 

    lumi_sf_23_PreBPix = 1.0  
    lumi_sf_23_PostBPix = 1.0
    lumi_sf_23 = 1.0
    lumi_sf_22 = 1.0 
    
    lumi_2022 = 34.6644
    lumi_2023_PreBPix = 17.9642 
    lumi_2023_PostBPix = 9.67674 
    lumi_2023 = lumi_2023_PreBPix + lumi_2023_PostBPix 
    lumi_total = lumi_2022 + lumi_2023_PreBPix + lumi_2023_PostBPix
 
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
    tree_sig.SetBranchStatus("event_HLT_HT_weight", 1) 
    tree_sig.SetBranchStatus("LLP*", 1) 
    tree_sig.SetBranchStatus("L1_prescale_weight", 1) 
    tree_sig.SetBranchStatus("HLT_prescale_weight", 1) 
    tree_sig.SetBranchStatus("puWeight*", 1) 
    tree_sig.SetBranchStatus("jet0_jet1_dPhi", 1)
    tree_sig.SetBranchStatus("Flag_METFilters_2022_2023_PromptReco", 1)
    tree_sig.SetBranchStatus("Pass_HLTDisplacedJet", 1) 
    tree_sig.SetBranchStatus("jet0_Tagged*", 1) 
    tree_sig.SetBranchStatus("jet1_Tagged*", 1) 
    
    train_frac = 0.2
    lj_train_cut = "(int(jet0_Pt * 1000) % 10) >= 8"
    sj_train_cut = "(int(jet1_Pt * 1000) % 10) >= 8"
 
    deltaPhi_cut = "(abs(jet0_jet1_dPhi) > 0.2)"

    # Lifetime Reweight

    reweight_llp0 = "pow ( {0} / {1}, 1 ) * exp( -LLP0_DecayCtau * 10. * ( 1.0/{2} - 1.0/{3} ) )".format(ctau_sample, ctau_target, ctau_target, ctau_sample)
    reweight_llp1 = "pow ( {0} / {1}, 1 ) * exp( -LLP1_DecayCtau * 10. * ( 1.0/{2} - 1.0/{3} ) )".format(ctau_sample, ctau_target, ctau_target, ctau_sample)
    
    # testing
#    if pileupweight_ == "nominal": reweight = "( puWeight * L1_prescale_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)
#    elif pileupweight_ == "pileupWeightUp": reweight = "( puWeightUp * L1_prescale_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)
#    elif pileupweight_ == "pileupWeightDown": reweight = "( puWeightDown * L1_prescale_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)

    # final
    if pileupweight_ == "nominal": reweight = "( puWeight * HLT_prescale_weight * L1_prescale_weight * event_HLT_HT_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)
    elif pileupweight_ == "pileupWeightUp": reweight = "( puWeightUp * HLT_prescale_weight * L1_prescale_weight * event_HLT_HT_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)
    elif pileupweight_ == "pileupWeightDown": reweight = "( puWeightDown * HLT_prescale_weight * L1_prescale_weight * event_HLT_HT_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)

    # no trigger
#    if pileupweight_ == "nominal": reweight = "( puWeight * HLT_prescale_weight * L1_prescale_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)
#    elif pileupweight_ == "pileupWeightUp": reweight = "( puWeightUp * HLT_prescale_weight * L1_prescale_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)
#    elif pileupweight_ == "pileupWeightDown": reweight = "( puWeightDown * HLT_prescale_weight * L1_prescale_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)

    # no pileup
    #reweight = "( HLT_prescale_weight * L1_prescale_weight * event_HLT_HT_weight * weight * {0} * {1})".format(reweight_llp0, reweight_llp1)

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

    hlt_pass =  "Pass_HLTDisplacedJet == 1 &&" if require_hlt_decision else ""
    tree_sig.Draw(
        "0.5 >> hist_sig_ljdc_"+ctau_target+"_"+sys+"_"+pileupweight_,
        " ({0}) * ( {6} Pass_PreSel == 1 && {1} && jet0_DepthTagCand == 1 && jet1_InclTagCand == 1 && jet0_scores_depth_LLPanywhere > {2} && jet1_scores_inc_train80 > {3} && {4} && {5})".format(
       reweight, lj_train_cut, lj_depth, lj_inc, deltaPhi_cut, trig_cut, hlt_pass
       )
    )

    tree_sig.Draw(
        "0.5 >> hist_sig_sjdc_"+ctau_target+"_"+sys+"_"+pileupweight_,
        " ({0}) * ( {6} Pass_PreSel == 1 && {1} && jet1_DepthTagCand == 1 && jet0_InclTagCand == 1 && jet1_scores_depth_LLPanywhere > {2} && jet0_scores_inc_train80 > {3} && {4} && {5})".format(
        reweight, sj_train_cut, sj_depth, sj_inc, deltaPhi_cut, trig_cut, hlt_pass
        )
    )

    nevents_sig_ljdc_temp = sigeff_sf_23 * hist_sig_ljdc.Integral() * SF_temp * 100. / (train_frac * lumi_sf_23)# 100 to convert from minituple % --> net fraction 
    nevents_sig_sjdc_temp = sigeff_sf_23 * hist_sig_sjdc.Integral() * SF_temp * 100. / (train_frac * lumi_sf_23)# 100 to convert from minituple % --> net fraction

    infile_sig.Close()

    return {
    "SIGLJDC_22": nevents_sig_ljdc_temp * lumi_2022 / lumi_total,
    "SIGLJDC_23PreBPix": nevents_sig_ljdc_temp * lumi_2023_PreBPix / lumi_total,
    "SIGLJDC_23PostBPix": nevents_sig_ljdc_temp * lumi_2023_PostBPix / lumi_total,
    "SIGSJDC_22": nevents_sig_sjdc_temp * lumi_2022 / lumi_total,
    "SIGSJDC_23PreBPix": nevents_sig_sjdc_temp * lumi_2023_PreBPix / lumi_total,
    "SIGSJDC_23PostBPix": nevents_sig_sjdc_temp * lumi_2023_PostBPix / lumi_total,
    }

# ------------------------------------------------------------------------------
def get_trigger_systematics(infilepath, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc, debug=False):
    MAX_TRIG = 2.0
    MAX_LEG_PCT = 50.0   # legs exceeding this are exclude from sys for now

    legs = [
        ("L1_lead",    "L1",    0, "pT", "jet0"),
        ("L1_sublead", "L1",    0, "pT", "jet1"),
        ("HLT1a_pT",  "HLT1a", 0, "pT", None),
        ("HLT1a_PTrk","HLT1a", 1, "PTrk", None),
        ("HLT1b_pT",  "HLT1b", 0, "pT", None),
        ("HLT1b_PTrk","HLT1b", 1, "PTrk", None),
        ("HLT2_pT",   "HLT2",  0, "pT", None),
        ("HLT2_PTrk", "HLT2",  1, "PTrk", None),
        ("HLT3a_pT",  "HLT3a", 0, "pT", None),
        ("HLT3a_PTrk","HLT3a", 1, "PTrk", None),
        ("HLT3a_DTrk","HLT3a", 2, "DTrk", None),
        ("HLT3b_pT",  "HLT3b", 0, "pT", None),
        ("HLT3b_PTrk","HLT3b", 1, "PTrk", None),
        ("HLT3b_DTrk","HLT3b", 2, "DTrk", None),
    ]
    regions = ["SIGLJDC_23PostBPix", "SIGSJDC_23PostBPix","SIGLJDC_23PreBPix", "SIGSJDC_23PreBPix", "SIGLJDC_22", "SIGSJDC_22"]

    per_filter_pct = {}   # raw, uncapped, unfiltered used for the printed table

    for label, fam, leg_idx, leg_name, single_jet in legs:

        if fam == "L1":
            nom_cut = "{0}_Tagged_{1}[{2}] == 1".format(single_jet, fam, leg_idx)
            var_cut = "{0}_Tagged_Varied_{1}[{2}] == 1".format(single_jet, fam, leg_idx)
        else:
            nom_cut = "(jet0_Tagged_{0}[{1}] == 1 && jet1_Tagged_{0}[{1}] == 1)".format(fam, leg_idx)

            gate0 = " && ".join("jet0_Tagged_{0}[{1}] == 1".format(fam, k) for k in range(leg_idx)) or "1"
            gate1 = " && ".join("jet1_Tagged_{0}[{1}] == 1".format(fam, k) for k in range(leg_idx)) or "1"

            var_cut = "(({0}) && jet0_Tagged_Varied_{2}[{3}] == 1) && (({1}) && jet1_Tagged_Varied_{2}[{3}] == 1)".format(
                gate0, gate1, fam, leg_idx
            )

        nominal = get_signal_yield(infilepath, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc,
                                    sys="trig_{0}_nom".format(label), trig_cut=nom_cut, require_hlt_decision=False)
        varied  = get_signal_yield(infilepath, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc,
                                    sys="trig_{0}_var".format(label), trig_cut=var_cut, require_hlt_decision=False)

        for region in regions:
            nom = nominal[region]
            var = varied[region]
            pct = 0.0 if nom <= 0 or nom < 0.1 else (var - nom) / nom * 100.0
            per_filter_pct["{0}_{1}".format(label, region)] = pct

            if debug:
                print("leg={0} region={1} nom={2} var={3} pct={4:.2f}%".format(label, region, nom, var, pct))

    hlt_legs = [
        "HLT1a_pT", "HLT1a_PTrk",
        "HLT1b_pT", "HLT1b_PTrk",
        "HLT2_pT", "HLT2_PTrk",
        "HLT3a_pT", "HLT3a_PTrk", "HLT3a_DTrk",
        "HLT3b_pT", "HLT3b_PTrk", "HLT3b_DTrk",
    ]

    combined = {}
    for region in regions:
        l1_lead_pct = per_filter_pct["L1_lead_{0}".format(region)]
        l1_sublead_pct = per_filter_pct["L1_sublead_{0}".format(region)]

        quad_sum_pct = 0.0
        if abs(l1_lead_pct) <= MAX_LEG_PCT:
            quad_sum_pct += l1_lead_pct ** 2
        else:
            print("WARNING: leg L1_lead region {0} pct={1:.2f}% exceeds {2}%, EXCLUDED from combined sys".format(
                region, l1_lead_pct, MAX_LEG_PCT))

        if abs(l1_sublead_pct) <= MAX_LEG_PCT:
            quad_sum_pct += l1_sublead_pct ** 2
        else:
            print("WARNING: leg L1_sublead region {0} pct={1:.2f}% exceeds {2}%, EXCLUDED from combined sys".format(
                region, l1_sublead_pct, MAX_LEG_PCT))

        for label in hlt_legs:
            leg_pct = per_filter_pct["{0}_{1}".format(label, region)]
            if abs(leg_pct) <= MAX_LEG_PCT:
                quad_sum_pct += leg_pct ** 2
            else:
                print("WARNING: leg {0} region {1} pct={2:.2f}% exceeds {3}%, EXCLUDED from combined sys".format(
                    label, region, leg_pct, MAX_LEG_PCT))

        val = 1.0 + (quad_sum_pct ** 0.5) / 100.0
        val = min(max(val, 1.0), MAX_TRIG)
        if val >= MAX_TRIG:
            print("WARNING LARGE COMBINED TRIGGER SYS:", region, val)
        combined[region] = val

    return combined, per_filter_pct
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
            # FINAL sideband update
            # sigeff_sf_23 = to update from 0.72 
            "lj_depth": 0.985,
            "lj_inc": 0.995,
            "sj_depth": 0.995,
            "sj_inc": 0.975,
            "bkg": {
                "2022": {"lj": 10.23, "sj": 0.23},
                "2023PreBPix": {"lj": 6.57, "sj": 1.23},
                "2023PostBPix": {"lj": 3.68, "sj": 1.82},
            },
        },
#            #WP 1 full opt. with tighter side bands and pT 3 bins
#            # sigeff_sf_23 = 0.71 
#            "lj_depth": 0.985,
#            "lj_inc": 0.995,
#            "sj_depth": 0.995,
#            "sj_inc": 0.975,
#            "bkg": {
#                2022: {"lj": 10.83, "sj": 0.000001},
#                2023: {"lj": 11.79, "sj": 2.95},
#            },
#        },
#            #WP 1 full opt. with tighter side bands
#            # sigeff_sf_23 = 0.71 
#            "lj_depth": 0.985,
#            "lj_inc": 0.995,
#            "sj_depth": 0.995,
#            "sj_inc": 0.975,
#            "bkg": {
#                2022: {"lj": 10.83, "sj": 0.0000001},
#                2023: {"lj": 11.52, "sj": 2.97},
#            },
#        },
#
#            #WP 1 full opt. with 400+GeV fix mistag rate method
#            "lj_depth": 0.985,
#            "lj_inc": 0.995,
#            "sj_depth": 0.995,
#            "sj_inc": 0.975,
#            "bkg": {
#                2022: {"lj": 4.46, "sj": 0.53},
#                2023: {"lj": 7.57, "sj": 1.83},
#            },
#        },
#        {
#            #WP 2 full opt. 2023 postBPix cuts only used
#            # sigeff_sf_23 = 0.29
#            "lj_depth": 0.995,
#            "lj_inc": 0.995,
#            "sj_depth": 0.995,
#            "sj_inc": 0.985,
#            "bkg": {
#                2022: {"lj": 0.55, "sj": 0.53},
#                2023: {"lj": 1.33, "sj": 1.14},
#            },
#        },
#        {
#            #WP 3 sig. eff. matched 2023 postBPix cuts only used
#            # sigeff_sf_23 = 1.0
#            "lj_depth": 0.995,
#            "lj_inc": 0.995,
#            "sj_depth": 0.995,
#            "sj_inc": 0.995,
#            "bkg": {
#                2022: {"lj": 6.23, "sj": 3.24},
#                2023: {"lj": 8.16, "sj": 3.77},
#            },
#        },
#   
    ]
    print("Reading in data tree... (this may take a few minutes)")
   
    infile_data_23PreBPix = [
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Cv1_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Cv2_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Cv3_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Cv4_scores.root",
    ]
    infile_data_23PostBPix = [
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Dv1_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2023Dv2_scores.root",
    ]
    infile_data_22 = [
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2022Dv1_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2022Ev1_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2022Fv1_scores.root",
    "/eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v5.6/minituple_data_2022Gv1_scores.root",
    ]
    
    tree_data_23PreBPix = make_chain("NoSel", infile_data_23PreBPix) 
    tree_data_23PostBPix = make_chain("NoSel", infile_data_23PostBPix) 
    tree_data_22   = make_chain("NoSel", infile_data_22)
    
    # Copy tree but only copy these branches
    tree_data_23PreBPix.SetBranchStatus("*", 0) 
    tree_data_23PreBPix.SetBranchStatus("Pass_PreSel", 1) 
    tree_data_23PreBPix.SetBranchStatus("jet*_DepthTagCand", 1) 
    tree_data_23PreBPix.SetBranchStatus("jet*_InclTagCand", 1) 
    tree_data_23PreBPix.SetBranchStatus("jet*_scores*", 1) 
    tree_data_23PreBPix.SetBranchStatus("jet*_DeepCSV*", 1)
    
    tree_data_23PostBPix.SetBranchStatus("*", 0) 
    tree_data_23PostBPix.SetBranchStatus("Pass_PreSel", 1) 
    tree_data_23PostBPix.SetBranchStatus("jet*_DepthTagCand", 1) 
    tree_data_23PostBPix.SetBranchStatus("jet*_InclTagCand", 1) 
    tree_data_23PostBPix.SetBranchStatus("jet*_scores*", 1) 
    tree_data_23PostBPix.SetBranchStatus("jet*_DeepCSV*", 1)
 
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
        nevents_sig_ljdc_23PreBPix = []
        nevents_sig_sjdc_23PreBPix = []
        nevents_sig_ljdc_23PostBPix = []
        nevents_sig_sjdc_23PostBPix = []
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
        tree_data_skim_23PreBPix = tree_data_23PreBPix.CopyTree("Pass_PreSel == 1")
        tree_data_skim_23PostBPix = tree_data_23PostBPix.CopyTree("Pass_PreSel == 1")
        tree_data_skim_22 = tree_data_22.CopyTree("Pass_PreSel == 1")
      
        nevents_bkg_ljdc_srpred_23PreBPix, nevents_bkg_sjdc_srpred_23PreBPix  = row["bkg"]["2023PreBPix"]["lj"], row["bkg"]["2023PreBPix"]["sj"]
        nevents_bkg_ljdc_srpred_23PostBPix, nevents_bkg_sjdc_srpred_23PostBPix  = row["bkg"]["2023PostBPix"]["lj"], row["bkg"]["2023PostBPix"]["sj"]
        nevents_bkg_ljdc_srpred_22, nevents_bkg_sjdc_srpred_22  = row["bkg"]["2022"]["lj"], row["bkg"]["2022"]["sj"] 
     
        # ----- Loop over Signal Lifetimes ----- #
        print( "Getting Event Counts...")

        for ctau_target in lifetimes:

            print( "\nCTau Target:", ctau_target )

            nominal = get_signal_yield(args.input, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc)

            nevents_sig_ljdc_23PreBPix.append( nominal["SIGLJDC_23PreBPix"] / SF_temp )
            nevents_sig_sjdc_23PreBPix.append( nominal["SIGSJDC_23PreBPix"] / SF_temp )
            nevents_sig_ljdc_23PostBPix.append( nominal["SIGLJDC_23PostBPix"] / SF_temp )
            nevents_sig_sjdc_23PostBPix.append( nominal["SIGSJDC_23PostBPix"] / SF_temp )
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
                    if debug: print( "negative:", nom )
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

            legs = [
                ("L1_lead",    "L1",    0, "pT"),
                ("L1_sublead", "L1",    0, "pT"),
                ("HLT1a_pT",  "HLT1a", 0, "pT"),
                ("HLT1a_PTrk","HLT1a", 1, "PTrk"),
                ("HLT1b_pT",  "HLT1b", 0, "pT"),
                ("HLT1b_PTrk","HLT1b", 1, "PTrk"),
                ("HLT2_pT",   "HLT2",  0, "pT"),
                ("HLT2_PTrk", "HLT2",  1, "PTrk"),
                ("HLT3a_pT",  "HLT3a", 0, "pT"),
                ("HLT3a_PTrk","HLT3a", 1, "PTrk"),
                ("HLT3a_DTrk","HLT3a", 2, "DTrk"),
                ("HLT3b_pT",  "HLT3b", 0, "pT"),
                ("HLT3b_PTrk","HLT3b", 1, "PTrk"),
                ("HLT3b_DTrk","HLT3b", 2, "DTrk"),
            ]
            
            trig_combined, trig_per_filter_pct = get_trigger_systematics(args.input, ctau_sample, ctau_target, lj_depth, lj_inc, sj_depth, sj_inc, debug=debug)
            
            print("--- Trigger SF variation by leg (ctau={0}) for {1}---".format(ctau_target, filetag))
            print("{0:<14}{1:>12}{2:>12}{3:>12}{4:>12}".format("Leg", "LJDC_23PreBPix", "SJDC_23PreBPix","LJDC_23PostBPix", "SJDC_23PostBPix", "LJDC_22", "SJDC_22"))
            for label, fam, leg_idx, leg_name in legs:
                row = [trig_per_filter_pct["{0}_SIG{1}".format(label, r)] for r in ["LJDC_23PreBPix", "SJDC_23PreBPix", "LJDC_23PostBPix", "SJDC_23PostBPix", "LJDC_22", "SJDC_22"]]
                print("{0:<14}{1:>11.2f}%{2:>11.2f}%{3:>11.2f}%{4:>11.2f}%".format(label, *row))

            # Replace test in template datacard
            output_file = template_datacard.replace("TEMPLATE", unique_filetag + "__" + ctau_target )

            print("Nevents LJDC 23PostBPix:", nominal["SIGLJDC_23PostBPix"])
            print("Nevents SJDC 23PostBPix:", nominal["SIGSJDC_23PostBPix"])
            print("Nevents LJDC 23PreBPix:", nominal["SIGLJDC_23PreBPix"])
            print("Nevents SJDC 23PreBPix:", nominal["SIGSJDC_23PreBPix"])
            print("Nevents LJDC 22:", nominal["SIGLJDC_22"])
            print("Nevents SJDC 22:", nominal["SIGSJDC_22"])

            replacements = {
                "SIGLJDC_23PostBPix": f"{nominal['SIGLJDC_23PostBPix']:.6e}", 
                "SIGSJDC_23PostBPix": f"{nominal['SIGSJDC_23PostBPix']:.6e}",
                "SIGLJDC_23PreBPix": f"{nominal['SIGLJDC_23PreBPix']:.6e}", 
                "SIGSJDC_23PreBPix": f"{nominal['SIGSJDC_23PreBPix']:.6e}",
                "SIGLJDC_22": f"{nominal['SIGLJDC_22']:.6e}", 
                "SIGSJDC_22": f"{nominal['SIGSJDC_22']:.6e}",

                "JERLJDC_23PostBPix": f"{jer['SIGLJDC_23PostBPix']:.6e}",
                "JERSJDC_23PostBPix": f"{jer['SIGSJDC_23PostBPix']:.6e}",
                "JERLJDC_23PreBPix": f"{jer['SIGLJDC_23PreBPix']:.6e}",
                "JERSJDC_23PreBPix": f"{jer['SIGSJDC_23PreBPix']:.6e}",
                "JERLJDC_22": f"{jer['SIGLJDC_22']:.6e}",
                "JERSJDC_22": f"{jer['SIGSJDC_22']:.6e}",

                "PULJDC_23PostBPix": f"{pu['SIGLJDC_23PostBPix']:.6e}",
                "PUSJDC_23PostBPix": f"{pu['SIGSJDC_23PostBPix']:.6e}",
                "PULJDC_23PreBPix": f"{pu['SIGLJDC_23PreBPix']:.6e}",
                "PUSJDC_23PreBPix": f"{pu['SIGSJDC_23PreBPix']:.6e}",
                "PULJDC_22": f"{pu['SIGLJDC_22']:.6e}",
                "PUSJDC_22": f"{pu['SIGSJDC_22']:.6e}",
                "TRIGLJDC_23PostBPix": f"{trig_combined['SIGLJDC_23PostBPix']:.6e}", 
                "TRIGSJDC_23PostBPix": f"{trig_combined['SIGSJDC_23PostBPix']:.6e}",
                "TRIGLJDC_23PreBPix": f"{trig_combined['SIGLJDC_23PreBPix']:.6e}", 
                "TRIGSJDC_23PreBPix": f"{trig_combined['SIGSJDC_23PreBPix']:.6e}",
                "TRIGLJDC_22": f"{trig_combined['SIGLJDC_22']:.6e}",
                "TRIGSJDC_22": f"{trig_combined['SIGSJDC_22']:.6e}",

                "BKGLJDC_23PostBPix": f"{nevents_bkg_ljdc_srpred_23PostBPix:.6e}", 
                "BKGSJDC_23PostBPix": f"{nevents_bkg_sjdc_srpred_23PostBPix:.6e}",
                "BKGLJDC_23PreBPix": f"{nevents_bkg_ljdc_srpred_23PreBPix:.6e}", 
                "BKGSJDC_23PreBPix": f"{nevents_bkg_sjdc_srpred_23PreBPix:.6e}",
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
        data["nevents_sig_ljdc_23PostBPix"] = nevents_sig_ljdc_23PostBPix
        data["nevents_sig_sjdc_23PostBPix"] = nevents_sig_sjdc_23PostBPix
        data["nevents_sig_ljdc_23PreBPix"] = nevents_sig_ljdc_23PreBPix
        data["nevents_sig_sjdc_23PreBPix"] = nevents_sig_sjdc_23PreBPix
        data["nevents_sig_ljdc_22"] = nevents_sig_ljdc_22
        data["nevents_sig_sjdc_22"] = nevents_sig_sjdc_22
        data["nevents_bkg_ljdc_23PostBPix"] = nevents_bkg_ljdc_srpred_23PostBPix
        data["nevents_bkg_sjdc_23PostBPix"] = nevents_bkg_sjdc_srpred_23PostBPix
        data["nevents_bkg_ljdc_23PreBPix"] = nevents_bkg_ljdc_srpred_23PreBPix
        data["nevents_bkg_sjdc_23PreBPix"] = nevents_bkg_sjdc_srpred_23PreBPix
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
        print( "LJDC 23 PostBPix:  ", nevents_sig_ljdc_23PostBPix )
        print( "SJDC 23 PostBPix:  ", nevents_sig_sjdc_23PostBPix )
        print( "LJDC 23 PreBPix:  ", nevents_sig_ljdc_23PreBPix )
        print( "SJDC 23 PreBPix:  ", nevents_sig_sjdc_23PreBPix )
        print( "LJDC 22:  ", nevents_sig_ljdc_22 )
        print( "SJDC 22:  ", nevents_sig_sjdc_22 )
        print( "BKG LJDC 23 PostBPix:  ", nevents_bkg_ljdc_srpred_23PostBPix )
        print( "BKG SJDC 23 PostBPix:  ", nevents_bkg_sjdc_srpred_23PostBPix )
        print( "BKG LJDC 23 PreBPix:  ", nevents_bkg_ljdc_srpred_23PreBPix )
        print( "BKG SJDC 23 PreBPix:  ", nevents_bkg_sjdc_srpred_23PreBPix )
        print( "BKG LJDC 22:  ", nevents_bkg_ljdc_srpred_22 )
        print( "BKG SJDC 22:  ", nevents_bkg_sjdc_srpred_22 )
     
        print( "--------------------------------------" )
        print( "Json file written to:", outfile_path )

if __name__ == '__main__':
    main()
