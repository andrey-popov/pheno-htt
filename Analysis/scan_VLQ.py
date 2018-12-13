#!/usr/bin/env python

"""Computes significance for model with vector-like quarks.

Performs a scan over the (mH, mQ) plane for given CP state and computes
a difference in significance with respect to pure SM+S model.  Results
of the scan can be stored in an .npz file.
"""

import argparse
import math
import os

import numpy as np

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt

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
        
        # Naive k-factors.  Set to the same values as for 2HDM.
        self.k_res = 2.
        self.k_int = math.sqrt(2. * 2.)
    
    
    def xsec_res(self, sqrt_s, alpha_s):
        """Compute cross section for resonant gg -> S -> tt."""
        
        s = sqrt_s ** 2
        
        if s <= 4 * self.mt ** 2:
            return 0.
        
        width = self.width_tt(self.cp, self.mass, s, g=self.g_tt)
        
        if self.cp == 'H':
            beta_power = 3
        else:
            beta_power = 1
        
        a = 3 * (alpha_s * self.gF * self.mt) ** 2 / (8192 * math.pi ** 3)
        loop_ampl = self.g_tt * self.loop_ampl_fermion(self.cp, s, mf=self.mt) \
            + self.num_vlq * self.g_vlq * self.loop_ampl_fermion(self.cp, s, mf=self.mass_vlq)
        denom = (s - self.mass ** 2) ** 2 + (width * self.mass) ** 2
        
        xsec = 2 * a * s ** 2 * self.beta(s) ** beta_power * self.g_tt ** 2 \
            * abs(loop_ampl) ** 2 / denom
        return self.to_pb(xsec * self.k_res)
    
    
    def xsec_int(self, sqrt_s, alpha_s):
        """Compute cross section for interference in gg -> S -> tt."""
        
        s = sqrt_s ** 2
        
        if s <= 4 * self.mt ** 2:
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
        
        loop_ampl = self.g_tt * self.loop_ampl_fermion(self.cp, s, mf=self.mt) \
            + self.num_vlq * self.g_vlq * self.loop_ampl_fermion(self.cp, s, mf=self.mass_vlq)
        propagator = s - self.mass ** 2 + 1j * width * self.mass
        
        xsec = a * b * beta ** beta_power * self.g_tt * (loop_ampl / propagator).real
        return self.to_pb(xsec * self.k_int)


if __name__ == '__main__':
    
    arg_parser = argparse.ArgumentParser(epilog=__doc__)
    arg_parser.add_argument(
        '--cp', default='H',
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
        higgs_masses = list(np.arange(350, 501, 10)) + list(np.arange(525, 1001, 25))
        vlq_masses = list(np.arange(200, 451, 10)) + list(np.arange(475, 601, 25)) \
            + list(np.arange(650, 1001, 50))
        
        grid = statscan.Grid(higgs_masses, vlq_masses)
        calc = statscan.StatCalc(None, args.bkg, args.lumi * 1e3)
        
        
        # Compute significance in SM+S
        sms_significance = {}
        
        for mH in higgs_masses:
            
            # SM+S is reproduced in the VLQ model when VLQ coupling is
            # set to zero
            parton_xsec = XSecVLQ(args.cp, mH, mH, g_vlq=0.)
            reco_mtt = RecoMtt(parton_xsec, resolution=args.resolution)
            calc.update_signal(reco_mtt)
            
            sms_significance[mH] = calc.significance()
        
        
        # Now perform the scan with VLQ
        for i, mH, mass_vlq in grid:
            
            parton_xsec = XSecVLQ(args.cp, mH, mass_vlq)
            reco_mtt = RecoMtt(parton_xsec, resolution=args.resolution)
            calc.update_signal(reco_mtt)
            
            significance = calc.significance()
            
            print('\033[1;34mResults for m{} = {:g}, mVLQ = {:g}:'.format(args.cp, mH, mass_vlq))
            print('  Significance: {}\n  Significance in SM+S: {}\033[0m'.format(
                significance, sms_significance[mH]
            ))
            
            # Instead of the significance, store the difference w.r.t.
            # SM+S.  CLs values will not be used, so put a placeholder.
            grid.set(i, significance - sms_significance[mH], 0.)
        
        
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
    fig = plt.figure()
    fig.patch.set_alpha(0.)
    axes = fig.add_subplot(111)
    
    xx_c, yy_c = np.meshgrid(
        statscan.PlotScan._centre_bins(grid.x), statscan.PlotScan._centre_bins(grid.y)
    )
    image = axes.pcolormesh(
        xx_c, yy_c, grid.significance.T, cmap='viridis',
        vmin=min(0., grid.significance.min())
    )
    fig.colorbar(
        image, fraction=0.05, pad=0.02,
        label='Diff. in significance w.r.t. SM+${}$ [$\\sigma$]'.format(args.cp),
        spacing='proportional'
    )
    
    axes.set_xlim(grid.x[0], grid.x[-1])
    axes.set_ylim(grid.y[1], grid.y[-1])
    
    # Mark line mH = 2 * mQ and label it.  The label is rotated
    # appropriately and offset a bit from the line
    axes.plot(
        [grid.x[0], grid.x[-1]], [grid.x[0] / 2, grid.x[-1] / 2],
        color='white', ls='dashed', lw=0.8
    )
    mark_pos_mH = grid.x[0] + 0.75 * (grid.x[-1] - grid.x[0])
    mark_pos_axes = (axes.transData + axes.transAxes.inverted()).transform(
        [mark_pos_mH, mark_pos_mH / 2]
    )
    mark_angle = axes.transData.transform_angles(
        np.array([math.degrees(0.5)]), np.array([[mark_pos_mH, mark_pos_mH / 2]])
    )[0]
    axes.text(
        mark_pos_axes[0], mark_pos_axes[1] + 0.01, '$m_{} = 2 m_Q$'.format(args.cp),
        size='small', color='white', transform=axes.transAxes,
        rotation=mark_angle, rotation_mode='anchor'
    )
    
    axes.set_xlabel('$m_{}$ [GeV]'.format(args.cp))
    axes.set_ylabel('$m_{Q}$ [GeV]')
    
    if args.lumi >= 1e3:
        lumi_text = '{:g} ab$^{{-1}}$'.format(args.lumi / 1e3)
    else:
        lumi_text = '{:g} fb$^{{-1}}$'.format(args.lumi)
    
    axes.text(
        0., 1.005, 'VLQ, CP-{} Higgs boson'.format('odd' if args.cp == 'A' else 'even'),
        ha='left', va='bottom', transform=axes.transAxes
    )
    axes.text(
        1., 1.005, 'Resolution {:g}%, $L = ${}'.format(args.resolution * 100, lumi_text),
        ha='right', va='bottom', transform=axes.transAxes
    )
    
    fig.savefig(args.output)
