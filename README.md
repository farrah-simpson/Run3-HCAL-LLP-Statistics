# Run3-HCAL-LLP-Statistics

Repository for Run 3 HCAL LLP statistical analysis

## Set up

Instructions and examples can be found in the central [Combine Documentation](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/?h=cmssw_14_1_0_pre4#combine-v10-recommended-version).

First time setup:
```
# Combine v10 -- recommended version
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git -c advice.detachedHead=false clone --depth 1 --branch v10.4.1 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
scramv1 b clean; scramv1 b -j$(nproc --ignore=2) # always make a clean build, with n - 2 cores on the system

# Clone this repository
git clone git@github.com:kennedykiley/Run3-HCAL-LLP-Statistics.git
cd Run3-HCAL-LLP-Statistics
git checkout -b <your branch name>
```

On every login:
```
cd CMSSW_14_1_0_pre4/src
cmsenv
cd HiggsAnalysis/CombinedLimit/Run3-HCAL-LLP-Statistics
```

## Run Combine Wrapper

The script `combine_wrapper.py` takes several inputs, including a template datacard (see `templates/`) and an input signal file. 

List input options in the combine wrapper:
```
python3 combine_wrapper.py -h
```
Example: This takes an input signal file and runs over multiple ctau target points, which are hardcoded in `combine_wrapper.py`. You must specify a filetag and signal point ctau (it is absolutely critical that the source ctau be correct). Optional arguments are the inclusive and depth score cuts that define the signal region. 
```
python3 combine_wrapper.py -t templates/v1/datacard_TEMPLATE.txt -i /eos/cms/store/group/phys_exotica/HCAL_LLP/MiniTuples/v4.1/minituple_HToSSTo4B_125_50_CTau3000_scores.root --filetag HToSSTo4B_125_50 --ctau 3000 --incl-score 0.9 --depth-score 0.8
```
This also makes .root files that are currently not used.

## Plot Limits

To plot limits, run:
```
mkdir plots/ # plots saved here

# Run limits for single mass point and multiple lifetimes:
python3 plot_limits.py <input_json_file> # makes typical brazil plot

# More for debugging: compare limits across multiple file configurations or signal mass points:
python3 plot_limits.py <output_filetag> <input_json_1> .... <input_json_N> # compares limits, signal yields, background yields
```
## How to Run Blinded Impact Plots

combineTool.py -M Impacts -d workspace.root -m 125 --doInitialFit -t -1 --rMin -100 --expectSignal 0 
combineTool.py -M Impacts -d workspace.root -m 125 --doFits --rMin -100 -t -1 --expectSignal 0 
combineTool.py -M Impacts -d  workspace.root-m 125 -o impacts.json
plotImpacts.py -i impact.json -o impacts

## To do

Improvements needed:
* More accurate background estimation (from Gillian)
* Application of correct SFs
* Integration of systematics beyond flat systematics in datacard
* More "finalized" SR cuts
* Beam halo rejection and VR outlier cut on btag score

Additional functionalities needed:
* Signal injection tests
* Asimov tests
* Pulls

