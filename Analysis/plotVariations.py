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
    
    figDir = 'figSyst'
    
    if not os.path.exists(figDir):
        os.makedirs(figDir)
    
    
    templatesFile = ROOT.TFile('ttbar.root')
    
    histNominal = templatesFile.Get('TT')
    nominal = hist_to_np(histNominal)
    
    # Extract binning from the nominal histogram
    binning = np.empty(histNominal.GetNbinsX() + 1)
    
    for bin in range(1, histNominal.GetNbinsX() + 2):
        binning[bin - 1] = histNominal.GetBinLowEdge(bin)
    
    
    # Plot individual variations, one per figure
    for systName, description, halfRange in [
        ('MttScale', 'Exp. scale in $m_{t\\bar t}$', 0.075),
        ('RenormScale', 'ME $\\mu_\\mathrm{R}$', 0.075),
        ('FactorScale', 'ME $\\mu_\\mathrm{F}$', 0.075),
        ('PDFAlphaS', '$\\alpha_s$ in PDF', 0.02),
        ('FSR', '$\\alpha_s$ in FSR', 0.075),
        ('MassT', '$m_t$', 0.075)
    ]:
        up = hist_to_np(templatesFile.Get('TT_{}Up'.format(systName))) / nominal - 1
        down = hist_to_np(templatesFile.Get('TT_{}Down'.format(systName))) / nominal - 1
        
        up *= 100
        down *= 100
        
        fig = plt.figure()
        axes = fig.add_subplot(111)
        
        axes.hist(binning[:-1], bins=binning, weights=up, histtype='step', label='Up')
        axes.hist(binning[:-1], bins=binning, weights=down, histtype='step', label='Down')
        axes.axhline(0., c='black', lw=0.8, ls='dashed')
        
        axes.margins(x=0.)
        axes.set_ylim(-halfRange * 100, halfRange * 100)
        
        axes.legend(loc='upper right')
        axes.set_xlabel('$m_{t\\bar t}$ [GeV]')
        axes.set_ylabel('Relative deviation from nominal [%]')
        
        axes.text(
            0.5, 1.005, description, size='large',
            ha='center', va='bottom', transform=axes.transAxes
        )
        
        fig.savefig(os.path.join(figDir, systName + '.pdf'))
        plt.close(fig)
    
    
    # All PDF variations are plotted in the same figure.  At the same
    # time identify small variations.
    smallPDFVars = []
    
    fig = plt.figure()
    axes = fig.add_subplot(111)
    
    nPDF = 30
    colourmap = plt.get_cmap('jet')
    
    for iPDF in range(nPDF):
        up = hist_to_np(templatesFile.Get('TT_PDF{}Up'.format(iPDF + 1))) / nominal - 1
        up *= 100
        axes.hist(
            binning[:-1], bins=binning, weights=up,
            color=colourmap(iPDF / (nPDF - 1)), histtype='step'
        )
        
        if np.max(up) < 0.1:
            smallPDFVars.append(iPDF + 1)
    
    axes.axhline(0., c='black', lw=0.8, ls='dashed')
    
    axes.set_ylim(-2, 2)
    axes.margins(x=0.)
    
    axes.set_xlabel('$m_{t\\bar t}$ [GeV]')
    axes.set_ylabel('Relative deviation from nominal [%]')
    
    axes.text(
        0.5, 1.005, 'PDF', size='large',
        ha='center', va='bottom', transform=axes.transAxes
    )
    
    fig.savefig(os.path.join(figDir, 'PDF.pdf'))
    plt.close(fig)
    
    print('PDF variations smaller than 0.1% everywhere:', smallPDFVars)
    
    templatesFile.Close()
