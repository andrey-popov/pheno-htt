#!/usr/bin/env python

"""Computes significance and CLs for SM+S and plots them.

Performs a scan over the (m, g) plane for each CP state.  Results of the
scan can be stored in an .npz file.
"""

import argparse
import os

import numpy as np

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from twohdm import XSecTwoHDM
from spectrum import RecoMtt
import statscan


if __name__ == '__main__':
    
    arg_parser = argparse.ArgumentParser(epilog=__doc__)
    arg_parser.add_argument(
        '--cp', default='A',
        help='Desired CP state'
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
    
    if args.cp not in {'A', 'H'}:
        raise RuntimeError('Cannot recognize CP state "{}".'.format(args.cp))
    
    if (args.bkg is None) == (args.from_file is None):
        raise RuntimeError('One and only one of options --bkg and --from-file must be given.')
    
    
    fig_dir = os.path.dirname(args.output)
    
    if fig_dir:
        try:
            os.makedirs(fig_dir)
        except FileExistsError:
            pass
    
    
    if not args.from_file:
        
        # Perform the scan if not reading results from a file
        mass_values = list(np.arange(350, 701, 25)) + list(np.arange(750, 1001, 50))
        
        if args.resolution < 0.15:
            if args.lumi < 200.:
                g_max = 3.
            elif args.lumi < 500.:
                g_max = 2.
            else:
                g_max = 1.5
        else:
            if args.lumi < 200.:
                g_max = 5.
            elif args.lumi < 500.:
                g_max = 4.
            else:
                g_max = 3.
        
        g_values = np.linspace(0.1, g_max, num=20)
        
        grid = statscan.Grid(mass_values, g_values)
        calc = statscan.StatCalc(None, args.bkg, args.lumi * 1e3)
        
        for i, mass, g in grid:
            
            width = XSecTwoHDM.width_tt(args.cp, mass, g=g)
            parton_xsec = XSecTwoHDM(mA=mass, wA=width, gA=g, mH=mass, wH=width, gH=g)
            
            if args.cp == 'A':
                parton_xsec.gH = 0.
            else:
                parton_xsec.gA = 0.
            
            reco_mtt = RecoMtt(parton_xsec, resolution=args.resolution)
            calc.update_signal(reco_mtt)
            
            significance = calc.significance()
            cls = calc.cls()
            
            print('\033[1;34mResults for {}, m = {:g}, g = {:g}:'.format(args.cp, mass, g))
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
    
    axes.set_xlabel('$m_{}$ [GeV]'.format(args.cp))
    axes.set_ylabel('$\\hat{g}_{\\Phi tt}$')
    
    axes.text(
        0., 1.005, 'SM+${}$'.format(args.cp),
        ha='left', va='bottom', transform=axes.transAxes
    )
    
    if args.lumi >= 1e3:
        lumi_text = '{:g} ab$^{{-1}}$'.format(args.lumi / 1e3)
    else:
        lumi_text = '{:g} fb$^{{-1}}$'.format(args.lumi)
    
    axes.text(
        1., 1.005, 'Resolution {:g}%, $L = ${}'.format(args.resolution * 100, lumi_text),
        ha='right', va='bottom', transform=axes.transAxes
    )
    
    fig.savefig(args.output)
