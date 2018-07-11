#!/usr/bin/env python

"""Computes significance and CLs for model with vector-like quarks.

Performs a scan over the (mH, mVLQ) plane for given CP state.  Results
of the scan can be stored in an .npz file.
"""

import math
import argparse
import os

import numpy as np

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from spectrum import PartonXSec, RecoMtt
import statscan


class XSecVLQ(PartonXSec):
    """Cross sections for gg -> S -> tt with vector-like quarks.
    
    A heavy Higgs boson of a fixed CP state is produced via gluon
    fusuon.  The loop is populated by top quarks and vector-like quarks
    (VLQ).  Multiple VLQ degenerate in mass are allowed.
    """
    
    def __init__(self, cp, mass, mass_vlq, g_tt=1., g_vlq=1., num_vlq=1):
        """Initialize from properties of the Higgs boson and VLQ.
        
        Arguments:
            cp:  CP state, 'A' or 'H'.
            mass:  Mass of the heavy Higgs boson, in GeV.
            mass_vlq:  Mass of vector-like quarks, in GeV.
            g_tt:  Reduced coupling of the heavy Higgs boson to top
                quarks.
            g_vlq:  Coupling of the heavy Higgs boson to vector-like
                quarks, normalized by the value of the coupling to top
                quarks in the SM.
            num_vlq:  Number of vector-like quarks with degenerate mass.
        """
        
        self._check_cp(cp)
        
        self.cp = cp
        self.mass = mass
        self.g_tt = g_tt
        self.g_vlq = g_vlq
        self.mass_vlq = mass_vlq
        self.num_vlq = num_vlq
        
        self.var_scale = self.width_tt(cp, mass, g=g_tt)
        
        # Naive k-factors
        self.k_res = 2.
        self.k_int = math.sqrt(2. * 1.6)
    
    
    def xsec_res(self, sqrt_s, alpha_s):
        """Compute cross section for resonant gg -> S -> tt."""
        
        s = sqrt_s ** 2
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        width = self.width_tt(self.cp, self.mass, s, g=self.g_tt)
        
        if self.cp == 'H':
            beta_power = 3
        else:
            beta_power = 1
        
        a = 3 * (alpha_s * self.gF * self.mt) ** 2 / (8192 * math.pi ** 3)
        loop_ampl = self.g_tt * self.loop_ampl(self.cp, s, mf=self.mt) \
            + self.num_vlq * self.g_vlq * self.loop_ampl(self.cp, s, mf=self.mass_vlq)
        denom = (s - self.mass ** 2) ** 2 + (width * self.mass) ** 2
        
        xsec = 2 * a * s ** 2 * self.beta(s) ** beta_power * self.g_tt ** 2 \
            * abs(loop_ampl) ** 2 / denom
        return self.to_pb(xsec * self.k_res)
    
    
    def xsec_int(self, sqrt_s, alpha_s):
        """Compute cross section for interference in gg -> S -> tt."""
        
        s = sqrt_s ** 2
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = self.beta(s)
        width = self.width_tt(self.cp, self.mass, s, g=self.g_tt)
        
        if self.cp == 'H':
            beta_power = 3
        else:
            beta_power = 1
        
        a = -alpha_s ** 2 * self.gF * self.mt ** 2 / (64 * math.pi * math.sqrt(2))
        
        # Factor dependent on z has been integrated
        b = math.log((1 + beta) / (1 - beta)) / beta
        
        loop_ampl = self.g_tt * self.loop_ampl(self.cp, s, mf=self.mt) \
            + self.num_vlq * self.g_vlq * self.loop_ampl(self.cp, s, mf=self.mass_vlq)
        propagator = s - self.mass ** 2 + 1j * width * self.mass
        
        xsec = a * b * beta ** beta_power * self.g_tt * (loop_ampl / propagator).real
        return self.to_pb(xsec * self.k_int)



if __name__ == '__main__':
    
    argParser = argparse.ArgumentParser(epilog=__doc__)
    argParser.add_argument(
        '--cp', default='A',
        help='Desired CP state'
    )
    argParser.add_argument(
        '-b', '--bkg', default=None,
        help='ROOT file with templates for SM tt'
    )
    argParser.add_argument(
        '-l', '--lumi', type=float, default=150.,
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
    
    if args.cp not in {'A', 'H'}:
        raise RuntimeError('Cannot recognize CP state "{}".'.format(args.cp))
    
    if (args.bkg is None) == (args.from_file is None):
        raise RuntimeError('One and only one of options --bkg and --from-file must be given.')
    
    
    fig_dir = os.path.dirname(args.output)
    
    if fig_dir and not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    
    
    if not args.from_file:
        
        # Perform the scan if not reading results from a file
        higgs_masses = np.arange(350, 1001, 50)
        vlq_masses = np.arange(1000, 2001, 100)
        
        grid = statscan.Grid(higgs_masses, vlq_masses)
        calc = statscan.StatCalc(None, args.bkg, args.lumi * 1e3)
        
        for i, mH, mass_vlq in grid:
            
            parton_xsec = XSecVLQ(args.cp, mH, mass_vlq)
            reco_mtt = RecoMtt(parton_xsec, resolution=args.resolution)
            calc.update_signal(reco_mtt)
            
            significance = calc.significance()
            cls = calc.cls()
            
            print('\033[1;34mResults for mH = {:g}, mVLQ = {:g}:'.format(mH, mass_vlq))
            print('  Significance: {}\n  CLs: {}\033[0m'.format(significance, cls))
            grid.set(i, significance, cls)
        
        
        # Save results of the scan if requested
        if args.save:
            scanner.save(args.save)
    
    else:
        # If reading scan results from a file, just load them
        grid = statscan.Grid.fromfile(args.from_file)
    
    
    # Plot results of the scan
    plotter = statscan.PlotScan(grid)
    fig, axes = plotter.draw()
    
    axes.set_xlabel('$m_{}$ [GeV]'.format(args.cp))
    axes.set_ylabel('$m_\\mathrm{VLQ}$ [GeV]')
    
    axes.text(
        0., 1.005, 'Vector-like quarks, CP-{}'.format('odd' if args.cp == 'A' else 'H'),
        ha='left', va='bottom', transform=axes.transAxex
    )
    
    if args.lumi >= 1e3:
        lumiText = '{:g} ab$^{{-1}}$'.format(args.lumi / 1e3)
    else:
        lumiText = '{:g} fb$^{{-1}}$'.format(args.lumi)
    
    axes.text(
        1., 1.005, 'Resolution {:g}%, $L = ${}'.format(args.resolution * 100, lumiText),
        ha='right', va='bottom', transform=axes.transAxes
    )
    
    fig.savefig(args.output)
