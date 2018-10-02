#!/usr/bin/env python

"""Computes significance and CLs for 2HDM and plots them.

Performs a scan over the (mA, tan(beta)) plane.  Results of the scan
can be stored in an .npz file.  Operates with benchmarks parameterized
with mA and tan(beta).
"""

import argparse
import math
import os

import numpy as np

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import twohdm
from spectrum import RecoMtt
import statscan


class XSecTwoHDMBenchmark(twohdm.XSecTwoHDM):
    """A class to compute cross section for gg -> S -> tt in 2HDM.
    
    Uses a fully specified benchmark defined with a ROOT file with
    parameters and parameterized with mA and tan(beta).
    """
    
    def __init__(self, mA, tanbeta, paramfile):
        
        paramfile = ROOT.TFile(paramfile)
        
        # Compute dependent parameters
        mH = paramfile.Get('m_H').Interpolate(mA, tanbeta)
        wA = paramfile.Get('width_A').Interpolate(mA, tanbeta)
        wH = paramfile.Get('width_H').Interpolate(mA, tanbeta)
        
        # Couplings in the decoupling limit are identical
        gA = gH = 1 / tanbeta
        
        
        # Set parameters of 2HDM
        super().__init__(mA, wA, gA, mH, wH, gH)
        self.tanbeta = tanbeta
        
        
        # Override scale over which the cross section changes
        self.var_scale = min(self.wA, self.wH)
        
        
        # Override dummy k-factors.  The k-factor for the SM tt
        # backround is set to 2 [1], and for the interference use the
        # geomentric mean of the k-factors for the background and the
        # resonant part
        # [1] https://github.com/andrey-popov/pheno-htt/issues/2
        # [2] Hespel et al., https://arxiv.org/abs/1606.04149
        k_bkg = 2.
        self.kA_res = paramfile.Get('kA_NNLO_13TeV').Interpolate(mA, tanbeta)
        self.kH_res = paramfile.Get('kH_NNLO_13TeV').Interpolate(mA, tanbeta)
        self.kA_int = math.sqrt(self.kA_res * k_bkg)
        self.kH_int = math.sqrt(self.kH_res * k_bkg)
        
        paramfile.Close()


if __name__ == '__main__':
    
    arg_parser = argparse.ArgumentParser(epilog=__doc__)
    arg_parser.add_argument(
        '-p', '--params', default=None,
        help='File with parameters that define the benchmark'
    )
    arg_parser.add_argument(
        '--label', default='', help='Label for the benchmark'
    )
    arg_parser.add_argument(
        '-b', '--bkg', default=None,
        help='ROOT file with templates for SM tt'
    )
    arg_parser.add_argument(
        '-l', '--lumi', type=float, default=150.,
        help='Target integrated luminosity, 1/fb'
    )
    arg_parser.add_argument(
        '-r', '--resolution', type=float, default=0.2,
        help='Relative resolution in mtt'
    )
    arg_parser.add_argument(
        '--save', default=None,
        help='Name of .npz file to store numeric results of the scan'
    )
    arg_parser.add_argument(
        '--from-file', dest='from_file', default=None,
        help='Name of .npz file with results of a scan'
    )
    arg_parser.add_argument(
        '-o', '--output', default='fig/significance.pdf',
        help='Name for output figure file'
    )
    args = arg_parser.parse_args()
    
    if (args.bkg is None) == (args.from_file is None):
        raise RuntimeError('One and only one of options --bkg and --from-file must be given.')
    
    if (args.params is None) == (args.from_file is None):
        raise RuntimeError('One and only one of options --params and --from-file must be given.')
    
    
    fig_dir = os.path.dirname(args.output)
    
    if fig_dir:
        try:
            os.makedirs(fig_dir)
        except FileExistsError:
            pass
    
    
    if not args.from_file:
        
        # Perform the scan if not reading results from a file
        mA_values = np.arange(350, 1001, 25)
        
        if args.resolution < 0.15:
            if args.lumi < 200.:
                max_tanbeta = 5.
            elif args.lumi < 500.:
                max_tanbeta = 6.
            else:
                max_tanbeta = 7.
        else:
            if args.lumi < 200.:
                max_tanbeta = 3.
            elif args.lumi < 500.:
                max_tanbeta = 4.
            else:
                max_tanbeta = 5.
        
        tanbeta_values = np.arange(0.5, max_tanbeta + 1e-3, 0.25)
        
        grid = statscan.Grid(mA_values, tanbeta_values)
        calc = statscan.StatCalc(None, args.bkg, args.lumi * 1e3)
        
        for i, mA, tanbeta in grid:
            
            parton_xsec = XSecTwoHDMBenchmark(mA, tanbeta, args.params)
            reco_mtt = RecoMtt(parton_xsec, resolution=args.resolution)
            calc.update_signal(reco_mtt)
            
            significance = calc.significance()
            cls = calc.cls()
            
            print('\033[1;34mResults for mA = {:g}, tan(beta) = {:g}:'.format(mA, tanbeta))
            print('  Significance: {}\n  CLs: {}\033[0m'.format(significance, cls))
            grid.set(i, significance, cls)
        
        
        # Save results of the scan if requested
        if args.save:
            scan_dir = os.path.dirname(args.save)
            
            if scan_dir:
                try:
                    os.makedirs(scan_dir)
                except FileExistsError:
                    pass
            
            grid.save(args.save)
    
    else:
        # If reading scan results from a file, just load them
        grid = statscan.Grid.fromfile(args.from_file)
    
    
    # Plot results of the scan
    plotter = statscan.PlotScan(grid)
    fig, axes = plotter.draw()
    
    axes.set_xlabel('$m_A$ [GeV]')
    axes.set_ylabel('$\\tan\/\\beta$')
    
    axes.text(0., 1.005, args.label, ha='left', va='bottom', transform=axes.transAxes)
    
    if args.lumi >= 1e3:
        lumi_text = '{:g} ab$^{{-1}}$'.format(args.lumi / 1e3)
    else:
        lumi_text = '{:g} fb$^{{-1}}$'.format(args.lumi)
    
    axes.text(
        1., 1.005, 'Resolution {:g}%, $L = ${}'.format(args.resolution * 100, lumi_text),
        ha='right', va='bottom', transform=axes.transAxes
    )
    
    fig.savefig(args.output)
