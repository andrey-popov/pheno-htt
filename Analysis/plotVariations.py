#!/usr/bin/env python

"""Plots systematic variations in SM tt."""

import os

import numpy as np

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from smoothTemplates import hist_to_np


if __name__ == '__main__':
    
    mpl.rc('xtick', top=True, direction='in')
    mpl.rc('ytick', right=True, direction='in')
    mpl.rc('axes', labelsize='large')
    mpl.rc('axes.formatter', limits=[-2, 4], use_mathtext=True)
    
    fig_dir = 'figSyst'
    
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    
    
    templates_file = ROOT.TFile('ttbar.root')
    
    hist_nominal = templates_file.Get('TT')
    nominal = hist_to_np(hist_nominal)
    
    # Extract binning from the nominal histogram
    binning = np.empty(hist_nominal.GetNbinsX() + 1)
    
    for bin in range(1, hist_nominal.GetNbinsX() + 2):
        binning[bin - 1] = hist_nominal.GetBinLowEdge(bin)
    
    
    # Plot individual variations, one per figure
    for syst_name, description, half_range in [
        ('MttScale', 'Exp. scale in $m_{t\\bar t}$', 0.075),
        ('RenormScale', 'ME $\\mu_\\mathrm{R}$', 0.075),
        ('FactorScale', 'ME $\\mu_\\mathrm{F}$', 0.075),
        ('PDFAlphaS', '$\\alpha_s$ in PDF', 0.02),
        ('FSR', '$\\alpha_s$ in FSR', 0.075),
        ('MassT', '$m_t$', 0.075)
    ]:
        up = hist_to_np(templates_file.Get('TT_{}Up'.format(syst_name))) / nominal - 1
        down = hist_to_np(templates_file.Get('TT_{}Down'.format(syst_name))) / nominal - 1
        
        up *= 100
        down *= 100
        
        fig = plt.figure()
        axes = fig.add_subplot(111)
        
        axes.hist(binning[:-1], bins=binning, weights=up, histtype='step', label='Up')
        axes.hist(binning[:-1], bins=binning, weights=down, histtype='step', label='Down')
        axes.axhline(0., c='black', lw=0.8, ls='dashed')
        
        axes.margins(x=0.)
        axes.set_ylim(-half_range * 100, half_range * 100)
        
        axes.legend(loc='upper right')
        axes.set_xlabel('$m_{t\\bar t}$ [GeV]')
        axes.set_ylabel('Relative deviation from nominal [%]')
        
        axes.text(
            0.5, 1.005, description, size='large',
            ha='center', va='bottom', transform=axes.transAxes
        )
        
        fig.savefig(os.path.join(fig_dir, syst_name + '.pdf'))
        plt.close(fig)
    
    
    # All PDF variations are plotted in the same figure.  At the same
    # time identify small variations.
    small_pdf_vars = []
    
    fig = plt.figure()
    axes = fig.add_subplot(111)
    
    num_pdf = 30
    colourmap = plt.get_cmap('jet')
    
    for iPDF in range(num_pdf):
        up = hist_to_np(templates_file.Get('TT_PDF{}Up'.format(iPDF + 1))) / nominal - 1
        up *= 100
        axes.hist(
            binning[:-1], bins=binning, weights=up,
            color=colourmap(iPDF / (num_pdf - 1)), histtype='step'
        )
        
        if np.max(up) < 0.1:
            small_pdf_vars.append(iPDF + 1)
    
    axes.axhline(0., c='black', lw=0.8, ls='dashed')
    
    axes.set_ylim(-2, 2)
    axes.margins(x=0.)
    
    axes.set_xlabel('$m_{t\\bar t}$ [GeV]')
    axes.set_ylabel('Relative deviation from nominal [%]')
    
    axes.text(
        0.5, 1.005, 'PDF', size='large',
        ha='center', va='bottom', transform=axes.transAxes
    )
    
    fig.savefig(os.path.join(fig_dir, 'PDF.pdf'))
    plt.close(fig)
    
    print('PDF variations smaller than 0.1% everywhere:', small_pdf_vars)
    
    templates_file.Close()
