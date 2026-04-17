import sys
import json
import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec


# --- Style settings (CMS-like) ---
plt.rcParams.update({
    "font.family": "sans-serif",
#    "font.sans-serif": ["Helvetica"], #, "Helvetica"],
    "mathtext.default": "regular",
    "axes.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True
})

# -------------------------------------------------------------------------------------------------
def get_data(infile):

	data_in = {}

	with open(infile) as f:
		data_in = json.load(f)

	mask = (np.array( data_in["limits_exp"]["50.0"] ) > 0)

	data_out = {}

	data_out["ctaus"] = np.array( data_in["ctaus"] )[mask] * 1.0e-3
	data_out["exp_median"] = np.array( data_in["limits_exp"]["50.0"] )[mask] #* 0.01
	data_out["exp_1s_low"] = np.array( data_in["limits_exp"]["16.0"] )[mask] #* 0.01
	data_out["exp_1s_high"] = np.array( data_in["limits_exp"]["84.0"] )[mask] #* 0.01
	data_out["exp_2s_low"] = np.array( data_in["limits_exp"][" 2.5"] )[mask] #* 0.01
	data_out["exp_2s_high"] = np.array( data_in["limits_exp"]["97.5"] )[mask] #* 0.01

	data_out["nevents_sig_ljdc_23"] = np.array( data_in["nevents_sig_ljdc_23"] )[mask] #* 0.01
	data_out["nevents_sig_sjdc_23"] = np.array( data_in["nevents_sig_sjdc_23"] )[mask] #* 0.01
	data_out["nevents_sig_ljdc_22"] = np.array( data_in["nevents_sig_ljdc_22"] )[mask] #* 0.01
	data_out["nevents_sig_sjdc_22"] = np.array( data_in["nevents_sig_sjdc_22"] )[mask] #* 0.01

	data_out["nevents_sig_ljdc"] = data_out["nevents_sig_ljdc_22"] + data_out["nevents_sig_ljdc_23"]
	data_out["nevents_sig_sjdc"] = data_out["nevents_sig_sjdc_22"] + data_out["nevents_sig_sjdc_23"]

	bkg_ljdc_23  = data_in["nevents_bkg_ljdc_23"]
	bkg_sjdc_23  = data_in["nevents_bkg_sjdc_23"]
	bkg_ljdc_22  = data_in["nevents_bkg_ljdc_22"]
	bkg_sjdc_22  = data_in["nevents_bkg_sjdc_22"]
	bkg_ljdc_total = (bkg_ljdc_22 + bkg_ljdc_23)
	bkg_sjdc_total = (bkg_sjdc_22 + bkg_sjdc_23)

	data_out["nevents_bkg_ljdc"] = np.full(len(data_out["ctaus"]), bkg_ljdc_total)
	data_out["nevents_bkg_sjdc"] = np.full(len(data_out["ctaus"]), bkg_sjdc_total)

	return data_out	

# -------------------------------------------------------------------------------------------------
def plot_single_limit(infile):

	# Load
	data = get_data(infile)

	#print(ctaus)
	#print(exp_median)

	fig, ax = plt.subplots(figsize=(5,5))

	# --- 2σ (yellow) and 1σ (green) bands ---
	ax.plot(data["ctaus"], data["exp_median"], 'k--', lw=2, label='Expected')
	ax.fill_between(data["ctaus"], data["exp_2s_low"], data["exp_2s_high"],
	                color='gold', alpha=0.5, label=r'$\pm2\sigma$')
	ax.fill_between(data["ctaus"], data["exp_1s_low"], data["exp_1s_high"],
	                color='limegreen', alpha=0.8, label=r'$\pm1\sigma$')

	# --- Median expected (black dashed) and observed (solid) ---
	#ax.plot(masses, obs, 'k-', lw=2, label='Observed')

	# --- Axes labels, limits, log scale ---
	ax.set_xlabel("CTau [m]", fontsize=13)
	ax.xaxis.label.set_horizontalalignment('right')
	ax.xaxis.set_label_coords(1.0, ax.xaxis.get_label().get_position()[1]-0.065)


	ax.set_ylabel(r"95% CL upper limit on BR(H$\to$SS)", fontsize=13)
	ax.yaxis.label.set_verticalalignment('top')
	ax.yaxis.set_label_coords(ax.yaxis.get_label().get_position()[0]-0.12, 0.66)

	ax.set_yscale("log")
	ax.set_xscale("log")
	ax.set_ylim(0.0005, 1.0)
	ax.grid(True, which="both", ls="--", lw=0.5, alpha=0.6)

	# --- CMS label and luminosity text ---
	#"""
	ax.text(0.0, 1.0, r"CMS", transform=ax.transAxes,
	        fontsize=16, fontweight='bold', va='bottom')
	ax.text(0.13, 1.0, r"Internal", transform=ax.transAxes,
	        fontsize=14, style='italic', va='bottom')
	ax.text(1.0, 1.0, r"63 fb$^{-1}$ (2022+2023) (13.6 TeV)", transform=ax.transAxes,
	        fontsize=12, ha='right', va='bottom')
	#"""
#60.42 fb-1
	# --- Legend ---
	ax.legend(loc="upper right", frameon=False, fontsize=11)

	plt.tight_layout()
	plt.subplots_adjust(top=0.92) 

	outfile = os.path.join("plots", infile.replace(".json", ".pdf").split("/")[-1])
	plt.savefig(outfile)

# -------------------------------------------------------------------------------------------------

def plot_limit_comparison(infiles, labels=None, outfile="plots/limit_comparison.pdf"):
    """
    Compare expected limits for multiple cut configurations on one plot.

    Parameters
    ----------
    infiles : list[str]
        List of json files.
    labels : list[str] or None
        Legend labels for each file. If None, use file basenames.
    outfile : str
        Output pdf path.
    """

    if labels is None:
        labels = [os.path.basename(f).replace(".json", "") for f in infiles]

    fig, ax = plt.subplots(figsize=(6, 5))

    for infile, label in zip(infiles, labels):
        data = get_data(infile)

        ax.plot(
            data["ctaus"],
            data["exp_median"],
            lw=2,
            label=label,
        )

    ax.set_xlabel("CTau [m]", fontsize=13)
    ax.xaxis.label.set_horizontalalignment('right')
    ax.xaxis.set_label_coords(1.0, ax.xaxis.get_label().get_position()[1] - 0.065)

    ax.set_ylabel(r"95% CL upper limit on BR(H$\to$SS)", fontsize=13)
    ax.yaxis.label.set_verticalalignment('top')
    ax.yaxis.set_label_coords(ax.yaxis.get_label().get_position()[0] - 0.12, 0.66)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.0005, 1.0)
    ax.grid(True, which="both", ls="--", lw=0.5, alpha=0.6)

    ax.text(0.0, 1.0, r"CMS", transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='bottom')
    ax.text(0.13, 1.0, r"Internal", transform=ax.transAxes,
            fontsize=14, style='italic', va='bottom')
    ax.text(1.0, 1.0, r"63 fb$^{-1}$ (2022+2023) (13.6 TeV)", transform=ax.transAxes,
            fontsize=12, ha='right', va='bottom')

    ax.legend(loc="best", frameon=False, fontsize=10)

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(outfile)
    plt.close(fig)

# -------------------------------------------------------------------------------------------------

def plot_ratio(infile1, infile2):
	# Load
	data1 = get_data(infile1)
	data2 = get_data(infile2)

	#print(ctaus)
	#print(exp_median)

	fig, ax = plt.subplots(figsize=(5,5))
    # --- ratios ---
	ctaus = data1["ctaus"]
	ratio_med = data1["exp_median"] / data2["exp_median"]
	ratio_1l = data1["exp_1s_low"]  / data2["exp_1s_low"]
	ratio_1h = data1["exp_1s_high"] / data2["exp_1s_high"]
	ratio_2l = data1["exp_2s_low"]  / data2["exp_2s_low"]
	ratio_2h = data1["exp_2s_high"] / data2["exp_2s_high"]
	ax.plot(ctaus, ratio_med, 'k--', lw=2, label='Expected Ratio')
	ax.fill_between(ctaus, ratio_2l, ratio_2h,
	            color='gold', alpha=0.5, label=r'$\pm2\sigma$ ratio')
	ax.fill_between(ctaus, ratio_1l, ratio_1h,
	            color='limegreen', alpha=0.8, label=r'$\pm1\sigma$ ratio')
	
	# --- Median expected (black dashed) and observed (solid) ---
	#ax.plot(masses, obs, 'k-', lw=2, label='Observed')
	
	# --- Axes labels, limits, log scale ---
	ax.set_xlabel("CTau [m]", fontsize=13)
	ax.xaxis.label.set_horizontalalignment('right')
	ax.xaxis.set_label_coords(1.0, ax.xaxis.get_label().get_position()[1]-0.065)


	ax.set_ylabel(r"Ratio of 95% CL upper limit on BR(H$\to$SS)", fontsize=13)
	ax.yaxis.label.set_verticalalignment('top')
	ax.yaxis.set_label_coords(ax.yaxis.get_label().get_position()[0]-0.12, 0.66)

	ax.set_xscale("log")
	ax.set_ylim(0.05, 2.0)
	ax.axhline(1.0, color="black", lw=1)

	# --- CMS label and luminosity text ---
	#"""
	ax.text(0.0, 1.0, r"CMS", transform=ax.transAxes,
	        fontsize=16, fontweight='bold', va='bottom')
	ax.text(0.13, 1.0, r"Internal", transform=ax.transAxes,
	        fontsize=14, style='italic', va='bottom')
	ax.text(1.0, 1.0, r"63 fb$^{-1}$ (2022+2023) (13.6 TeV)", transform=ax.transAxes,
	        fontsize=12, ha='right', va='bottom')
	#"""

	# --- Legend ---
	ax.legend(loc="upper right", frameon=False, fontsize=11)

	plt.tight_layout()
	plt.subplots_adjust(top=0.92) 

	outfile = os.path.join("plots", infile1.replace(".json", "_ratio.pdf").split("/")[-1])
	plt.savefig(outfile)


def plot_multi_limit(infiles):

	# Load
	data = {}
	filetags = []
	for infile in infiles:
		filetag = infile.replace(".json", "").split("/")[-1]
		filetags.append(filetag)
		data[filetag] = get_data(infile)

	#print(ctaus)
	#print(exp_median)

	fig, ax = plt.subplots(figsize=(10,8))

	color_map_name = "summer_r"
	cmap = get_cmap(color_map_name, len(filetags)) # plasma, viridis

	i = -1
	for filetag in filetags:
		i += 1

		if len(data[filetag]["ctaus"]) < 5: continue
		print(filetag, data[filetag]["exp_median"])
		continue
		
		# --- 2σ (yellow) and 1σ (green) bands ---
		#ax.fill_between(data["ctaus"], data["exp_2s_low"], data["exp_2s_high"],
		#                color='gold', alpha=0.5, label=r'$\pm2\sigma$ expected')
		#ax.fill_between(data["ctaus"], data["exp_1s_low"], data["exp_1s_high"],
		#                color='limegreen', alpha=0.8, label=r'$\pm1\sigma$ expected')

		# --- Median expected (black dashed) and observed (solid) ---
		ax.plot(data[filetag]["ctaus"], data[filetag]["exp_median"], color=cmap(i), lw=2, label=filetag)
		#ax.plot(masses, obs, 'k-', lw=2, label='Observed')

	quit()
	# --- Axes labels, limits, log scale ---
	ax.set_xlabel("CTau [m]", fontsize=13)
	ax.set_ylabel("95% CL upper limit on BR(H->XX * X->bb)", fontsize=13)
	ax.set_yscale("log")
	ax.set_xscale("log")
	#ax.set_ylim(0.0001, 1000.0)
	ax.grid(True, which="both", ls="--", lw=0.5, alpha=0.6)

	# --- CMS label and luminosity text ---
	"""
	ax.text(0.03, 0.95, r"CMS", transform=ax.transAxes,
	        fontsize=16, fontweight='bold', va='top')
	ax.text(0.17, 0.95, r"Preliminary", transform=ax.transAxes,
	        fontsize=13, style='italic', va='top')
	ax.text(0.97, 0.95, r"XX fb$^{-1}$ (13.6 TeV)", transform=ax.transAxes,
	        fontsize=11, ha='right', va='top')
	"""

	# --- Legend ---
	ax.legend(loc="upper right", frameon=False, fontsize=11)

	plt.tight_layout()

	outfile = os.path.join("plots", "multi.png")
	plt.savefig(outfile)	

# -------------------------------------------------------------------------------------------------
def plot_multi_limit_debug(outfiletag, infiles):

	# Load
	data = {}
	filetags = []
	for infile in infiles:
		filetag = infile.replace(".json", "").split("/")[-1]
		filetags.append(filetag)
		data[filetag] = get_data(infile)

	#print(ctaus)
	#print(exp_median)

	fig = plt.figure(figsize=(10, 8))
	gs = GridSpec(2, 1, height_ratios=[2, 1], figure=fig, hspace=0.0)  # 2:1 ratio

	ax_top = fig.add_subplot(gs[0])
	ax_bottom_L = fig.add_subplot(gs[1], sharex=ax_top)
	ax_bottom_R = ax_bottom_L.twinx()

	color_map_name = "summer_r"
	cmap = get_cmap(color_map_name, len(filetags)) # plasma, viridis

	i = -1
	for filetag in filetags:
		i += 1

		#if len(data[filetag]["ctaus"]) < 5: continue

		
		# --- 2σ (yellow) and 1σ (green) bands ---
		#ax.fill_between(data["ctaus"], data["exp_2s_low"], data["exp_2s_high"],
		#                color='gold', alpha=0.5, label=r'$\pm2\sigma$ expected')
		#ax.fill_between(data["ctaus"], data["exp_1s_low"], data["exp_1s_high"],
		#                color='limegreen', alpha=0.8, label=r'$\pm1\sigma$ expected')

		# --- Median expected (black dashed) and observed (solid) ---
		ax_top.plot(data[filetag]["ctaus"], data[filetag]["exp_median"], color=cmap(i), lw=2, label=filetag)

		if i == 0:
			ax_bottom_L.plot([],[], color='k', linestyle="-.", lw=2, label='Signal LJDC')
			ax_bottom_L.plot([],[], color='k', linestyle=":", lw=2, label='Signal SJDC')
			ax_bottom_L.plot([],[], color='k', linestyle="-", lw=2, label='Background LJDC')
			ax_bottom_L.plot([],[], color='k', linestyle="--", lw=2, label='Background LJDC')

		ax_bottom_L.plot(data[filetag]["ctaus"], data[filetag]["nevents_sig_ljdc"], color=cmap(i), linestyle="-.", lw=2)# , label='LJDC')
		ax_bottom_L.plot(data[filetag]["ctaus"], data[filetag]["nevents_sig_sjdc"], color=cmap(i), linestyle=":", lw=2)#, label='SJDC')

		ax_bottom_R.plot(data[filetag]["ctaus"], data[filetag]["nevents_bkg_ljdc"], color=cmap(i), linestyle="-", lw=2)#, label='Data LJDC')
		ax_bottom_R.plot(data[filetag]["ctaus"], data[filetag]["nevents_bkg_sjdc"], color=cmap(i), linestyle="--", lw=2)#, label='Data SJDC')

		#break

		#ax.plot(masses, obs, 'k-', lw=2, label='Observed')

	# --- Axes labels, limits, log scale ---
	#ax_top.set_xlabel("CTau [m]", fontsize=13)
	ax_top.set_ylabel("95% CL upper limit on BR(H->XX * X->bb)", fontsize=13)
	ax_top.set_yscale("log")
	ax_top.set_xscale("log")
	#ax.set_ylim(0.0001, 1000.0)
	ax_top.grid(True, which="both", ls="--", lw=0.5, alpha=0.6)

	ax_bottom_L.set_xlabel("CTau [m]", fontsize=13)
	ax_bottom_L.set_ylabel("N Signal Events")
	ax_bottom_L.set_yscale("log")
	ax_bottom_R.set_ylabel("N Background Events")
	ax_bottom_R.set_yscale("log")

	# --- CMS label and luminosity text ---
	"""
	ax.text(0.03, 0.95, r"CMS", transform=ax.transAxes,
	        fontsize=16, fontweight='bold', va='top')
	ax.text(0.17, 0.95, r"Preliminary", transform=ax.transAxes,
	        fontsize=13, style='italic', va='top')
	ax.text(0.97, 0.95, r"XX fb$^{-1}$ (13.6 TeV)", transform=ax.transAxes,
	        fontsize=11, ha='right', va='top')
	"""

	# --- Legend ---
	ax_top.legend(loc="upper right", frameon=False, fontsize=11)
	ax_bottom_L.legend(loc="upper right", frameon=False, fontsize=11)
	#ax_bottom_R.legend(loc="center right", frameon=False, fontsize=11)

	plt.tight_layout()

	outfile = os.path.join("plots", outfiletag+".png")
	print("Saving figure to:", outfile )
	plt.savefig(outfile)	

# -------------------------------------------------------------------------------------------------
def main():

    #single limit plot
    #	if len(sys.argv) == 2: 
    #		infile = sys.argv[1]
    #		plot_single_limit(infile)
    #		return
    #	else: 
    #		filetag = sys.argv[1]
    #		infiles = sys.argv[2:]
    #		plot_multi_limit_debug(filetag, infiles)
    
    #limit ratio plot
    #    plot_ratio("output/HToSSTo4B_125_50_inc0.9_depth0.8.json","output/HToSSTo4B_125_50_inc0.9_depth0.8_statonly.json")
    
    #limit comparison plot

    filetag = "HToSSTo4B_125_50" 
    
    outdir = "output/"
    labels = [
    #    "LJDC 0.965/0.695, SJDC 0.965/0.315",
        "LJDC 0.965/0.845, SJDC 0.975/0.375",
        "LJDC 0.975/0.415, SJDC 0.975/0.415",
    ]
    
    
    infiles = [
    #    os.path.join(outdir, "{0}_inc0.695_0.315_depth0.965_0.965.json".format(filetag)),
        os.path.join(outdir, "{0}_inc0.845_0.375_depth0.965_0.975.json".format(filetag)),
        os.path.join(outdir, "{0}_inc0.415_0.415_depth0.975_0.975.json".format(filetag)),
    ]
    
    plot_limit_comparison(infiles, labels, outfile="plots/limit_comparison.pdf")

# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	main()

