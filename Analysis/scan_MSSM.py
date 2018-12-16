#!/usr/bin/env python

"""Computes significance and CLs for MSSM and plots them.

Performs a scan over the (mA, tan(beta)) plane.  Results of the scan
can be stored in an .npz file.
"""

import argparse
import itertools
import math
import os

import numpy as np

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from spectrum import RecoMtt, PartonXSec
import statscan


class XSecMSSM(PartonXSec):
    """Cross section for gg -> S -> tt in MSSM.
    
    Takes into account the contribution from two stop quarks in the
    ggH coupling.  They do not contribute to the ggA vertex, though.
    """
    
    def __init__(self, mA, tanbeta, paramfile):
        """Initialize from mA, tan(beta) and file with parameters."""
        
        paramfile = ROOT.TFile(paramfile)
        
        self.mA = mA
        self.mH = paramfile.Get('m_H').Interpolate(mA, tanbeta)
        
        self.wA = paramfile.Get('width_A').Interpolate(mA, tanbeta)
        self.wH = paramfile.Get('width_H').Interpolate(mA, tanbeta)
        
        # Couplings to top quarks
        self.gA_top = 1 / tanbeta
        alpha = paramfile.Get('alphaA').Interpolate(mA, tanbeta)
        self.gH_top = math.sin(alpha) / (tanbeta / math.sqrt(1 + tanbeta ** 2))
        
        # Masses of stop quarks
        self.m_stop1 = paramfile.Get('mstop1A').Interpolate(mA, tanbeta)
        self.m_stop2 = paramfile.Get('mstop2A').Interpolate(mA, tanbeta)
        
        # Couplings of the CP-even Higgs boson to stop quarks.  The
        # CP-odd state does not couple to them.
        self.gH_stop1 = paramfile.Get('gstop11H').Interpolate(mA, tanbeta)
        self.gH_stop2 = paramfile.Get('gstop22H').Interpolate(mA, tanbeta)
        
        # Set k-factors.  The k-factor for the SM tt backround is set to
        # 2 [1], and for the interference use the geomentric mean of the
        # k-factors for the background and the resonant part
        # [1] https://github.com/andrey-popov/pheno-htt/issues/2
        # [2] Hespel et al., https://arxiv.org/abs/1606.04149
        k_bkg = 2.
        self.kA_res = paramfile.Get('kA_NNLO_13TeV').Interpolate(mA, tanbeta)
        self.kH_res = paramfile.Get('kH_NNLO_13TeV').Interpolate(mA, tanbeta)
        self.kA_int = math.sqrt(self.kA_res * k_bkg)
        self.kH_int = math.sqrt(self.kH_res * k_bkg)
        
        paramfile.Close()
        
        
        # Set scale over which the cross section changes
        self.var_scale = min(self.wA, self.wH)
    
    
    @staticmethod
    def loop_ampl_scalar(cp, s, m):
        """Compute scalar loop amplitude.
        
        Arguments:
            cp:  CP state of the Higgs boson, 'A' or 'H'.
            s:  Mandelstam s variable, in GeV^2.
            m:  Mass of the scalar running in the loop, in GeV.
        
        Return value:
            Computed amplitude.  It is a complex number.
        """
        
        if cp != 'H':
            raise NotImplementedError
        
        tau = s / (2 * m) ** 2
        
        if tau > 1:
            beta = math.sqrt(1 - 1 / tau)
            f = -0.25 * (math.log((1 + beta) / (1 - beta)) - math.pi * 1j) ** 2
        else:
            f = math.asin(math.sqrt(tau)) ** 2
        
        return -(tau - f) / tau ** 2
    
    
    def xsec_res(self, sqrt_s, alpha_s):
        """Compute cross section for resonant gg -> S -> tt."""
        
        s = sqrt_s ** 2
        
        if s <= 4 * self.mt ** 2:
            return 0.
        
        prefactor = 3 * (alpha_s * self.gF * self.mt) ** 2 / (8192 * math.pi ** 3)
        
        sum_scalars = 0.
        beta = self.beta(s)
        
        # Contributions from CP-odd state
        denom = (s - self.mA ** 2) ** 2 + (self.wA * self.mA) ** 2
        ampl = self.gA_top ** 2 * self.loop_ampl_fermion('A', s)
        sum_scalars += self.kA_res * beta * abs(ampl) ** 2 / denom
        
        # Contribution from CP-even state
        denom = (s - self.mH ** 2) ** 2 + (self.wH * self.mH) ** 2
        ampl = self.gH_top ** 2 * self.loop_ampl_fermion('H', s)
        
        for m_stop, gH_stop in [(self.m_stop1, self.gH_stop1), (self.m_stop2, self.gH_stop2)]:
            ampl += self.gH_top * gH_stop * (self.mt / m_stop) ** 2 / 2 * \
                self.loop_ampl_scalar('H', s, m_stop)
        
        sum_scalars += self.kH_res * beta ** 3 * abs(ampl) ** 2 / denom
        
        return self.to_pb(2 * prefactor * s ** 2 * sum_scalars)
    
    
    def xsec_int(self, sqrt_s, alpha_s):
        """Compute cross section for interference in gg -> S -> tt."""
        
        s = sqrt_s ** 2
        
        if s <= 4 * self.mt ** 2:
            return 0.
        
        prefactor = - alpha_s ** 2 * self.gF * self.mt ** 2 / (64 * math.sqrt(2) * math.pi)
        
        # Integral of the factor depending on z
        beta = self.beta(s)
        integral = 2 / beta * math.atanh(beta)
        
        sum_scalars = 0.
        
        # Contribution from CP-odd state
        denom = s - self.mA ** 2 + 1j * self.wA * self.mA
        ampl = self.gA_top ** 2 * self.loop_ampl_fermion('A', s)
        sum_scalars += self.kA_int * beta * ampl / denom
        
        # Contribution from CP-even state
        denom = s - self.mH ** 2 + 1j * self.wH * self.mH
        ampl = self.gH_top ** 2 * self.loop_ampl_fermion('H', s)
        
        for m_stop, gH_stop in [(self.m_stop1, self.gH_stop1), (self.m_stop2, self.gH_stop2)]:
            ampl += self.gH_top * gH_stop * (self.mt / m_stop) ** 2 / 2 * \
                self.loop_ampl_scalar('H', s, m_stop)
        
        sum_scalars += self.kH_int * beta ** 3 * ampl / denom
        
        return self.to_pb(prefactor * integral * sum_scalars.real)


if __name__ == '__main__':
    
    arg_parser = argparse.ArgumentParser(epilog=__doc__)
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
        
        tanbeta_values = np.arange(0.75, max_tanbeta + 1e-3, 0.25)
        
        grid = statscan.Grid(mA_values, tanbeta_values)
        calc = statscan.StatCalc(None, args.bkg, args.lumi * 1e3)
        
        for i, mA, tanbeta in grid:
            
            parton_xsec = XSecMSSM(mA, tanbeta, 'params/MSSM.root')
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
    
    axes.text(0., 1.005, 'MSSM', ha='left', va='bottom', transform=axes.transAxes)
    
    if args.lumi >= 1e3:
        lumi_text = '{:g} ab$^{{-1}}$'.format(args.lumi / 1e3)
    else:
        lumi_text = '{:g} fb$^{{-1}}$'.format(args.lumi)
    
    axes.text(
        1., 1.005, 'Resolution {:g}%, $L = ${}'.format(args.resolution * 100, lumi_text),
        ha='right', va='bottom', transform=axes.transAxes
    )
    
    
    # Mark the region with large ggH k-factors
    axes.set_xlim(*axes.get_xlim())
    axes.set_ylim(*axes.get_ylim())
    
    paramfile = ROOT.TFile('params/MSSM.root')
    hist_kfactor = paramfile.Get('kH_NNLO_13TeV')
    hist_kfactor.SetDirectory(None)
    paramfile.Close()
    
    x, y = [], []
    
    for bin in range(1, hist_kfactor.GetNbinsX() + 2):
        x.append(hist_kfactor.GetXaxis().GetBinLowEdge(bin))
    
    for bin in range(1, hist_kfactor.GetNbinsY() + 2):
        y.append(hist_kfactor.GetYaxis().GetBinLowEdge(bin))
    
    xx, yy = np.meshgrid(x, y)
    kfactor = np.empty_like(xx)
    
    for ix, iy in itertools.product(
        range(hist_kfactor.GetNbinsX()), range(hist_kfactor.GetNbinsY())
    ):
        kfactor[iy, ix] = hist_kfactor.GetBinContent(ix + 1, iy + 1)
    
    axes.contourf(
        xx, yy, kfactor, [10., math.inf],
        colors='none', hatches=['///'], zorder=1.5
    )
    
    fig.savefig(args.output)
