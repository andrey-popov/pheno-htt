#!/usr/bin/env python

"""Computes significance and CLs for hMSSM and plots them.

Performs a scan over the (mA, tan(beta)) plane.  Results of the scan
can be stored in an .npz file.
"""

import argparse
import json
import os

import numpy as np

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import hmssm
from spectrum import RecoMtt
import statscan


if __name__ == '__main__':
    
    argParser = argparse.ArgumentParser(epilog=__doc__)
    argParser.add_argument(
        '--binning', default='binning.json',
        help='JSON file with binning for reconstructed mtt'
    )
    argParser.add_argument(
        '--bkg', default='ttbar.root',
        help='ROOT file with templates for SM tt'
    )
    argParser.add_argument(
        '--lumi', type=float, default=150.,
        help='Target integrated luminosity, 1/fb'
    )
    argParser.add_argument(
        '-r', '--resolution', type=float, default=0.2,
        help='Relative resolution in mtt'
    )
    argParser.add_argument(
        '--save', default=None,
        help='Name of .npz file to store numeric results of the scan'
    )
    argParser.add_argument(
        '--from-file', dest='from_file', default=None,
        help='Name of .npz file with results of a scan'
    )
    argParser.add_argument(
        '-o', '--output', default='fig/significance.pdf',
        help='Name for output figure file'
    )
    args = argParser.parse_args()
    
    
    figDir = os.path.dirname(args.output)
    
    if figDir and not os.path.exists(figDir):
        os.makedirs(figDir)
    
    
    if not args.from_file:
        
        # Perform the scan if not reading results from a file
        with open(args.binning) as f:
            binning = np.array(json.load(f), dtype=np.float64)
        
        bkgfile = ROOT.TFile(args.bkg)
        
        
        mA_values = np.arange(350, 1001, 25)
        tanbeta_values = np.arange(0.75, 5.1, 0.25)
        
        grid = statscan.Grid(mA_values, tanbeta_values)
        
        for i, mA, tanbeta in grid:
            
            parton_xsec = hmssm.XSecHMSSM(mA, tanbeta, 'hMSSM_13TeV.root')
            reco_mtt = RecoMtt(parton_xsec, resolution=args.resolution)
            calc = statscan.StatCalc(reco_mtt, binning, bkgfile, args.lumi * 1e3)
            
            significance = calc.significance()
            cls = calc.cls()
            
            print('\033[1;34mResults for mA = {:g}, tan(beta) = {:g}:'.format(mA, tanbeta))
            print('  Significance: {}\n  CLs: {}\033[0m'.format(significance, cls))
            grid.set(i, significance, cls)
        
        
        # Save results of the scan if requested
        if args.save:
            scanner.save(args.save)
        
        bkgfile.Close()
    
    else:
        # If reading scan results from a file, just load them
        grid = statscan.Grid.fromfile(args.from_file)
    
    
    # Plot results of the scan
    plotter = statscan.PlotScan(grid)
    fig, axes = plotter.draw()
    
    axes.set_xlabel('$m_A$ [GeV]')
    axes.set_ylabel('$\\tan\/\\beta$')
    
    axes.text(0., 1.005, 'hMSSM', ha='left', va='bottom', transform=axes.transAxes)
    
    if args.lumi >= 1e3:
        lumiText = '{:g} ab$^{{-1}}$'.format(args.lumi / 1e3)
    else:
        lumiText = '{:g} fb$^{{-1}}$'.format(args.lumi)
    
    axes.text(
        1., 1.005, 'Resolution {:g}%, $L = ${}'.format(args.resolution * 100, lumiText),
        ha='right', va='bottom', transform=axes.transAxes
    )
    
    fig.savefig(args.output)
