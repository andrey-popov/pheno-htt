#!/usr/bin/env python

"""Suppresses fluctuations in systematic variations."""

import os

import numpy as np

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


def hist_to_np(hist, full=False):
    """Convert ROOT histogram content into a NumPy array.
    
    Arguments:
        hist:  ROOT histogram.
        full:  Indicates whether only bin contents are extracted or also
            the binning and errors.
    
    Return value:
        A NumPy array of bin contents if full is false, otherwise a
        tuple of arrays with bin edges, contents, and errors.
    """
    
    num_bins = hist.GetNbinsX()
    content = np.empty(num_bins)
    
    for bin in range(1, num_bins + 1):
        content[bin - 1] = hist.GetBinContent(bin)
    
    if not full:
        return content
    else:
        binning = np.empty(num_bins + 1)
        errors = np.empty(num_bins)
        
        for bin in range(1, num_bins + 1):
            binning[bin - 1] = hist.GetBinLowEdge(bin)
            errors[bin - 1] = hist.GetBinError(bin)
        
        binning[-1] = hist.GetBinLowEdge(num_bins + 1)
        
        return binning, content, errors


class Smoother:
    """A class to smooth fluctuations in a systematic variations.
    
    Smooths up and down relative deviations using LOWESS algorithm [1].
    The relative deviations are assumed to be symmetric in shape but
    allowed to scale differently.
    [1] https://en.wikipedia.org/wiki/Local_regression
    """
    
    def __init__(self, nominal, up, down, weights):
        """Initializer from reference templates.
        
        Arguments:
            nominal, up, down:  Nominal template and templates that
                describe up and down variations.
            weights:  Relative weights for bins of the templates.
        """
        
        self.nominal = nominal
        self.up = up
        self.down = down
        self.weights = weights
        
        self.smooth_averaged_deviation = None
    
    
    @staticmethod
    def lowess(y, external_weights, bandwidth):
        """Peform LOWESS smoothing.
        
        Arguments:
            y:  Input sequence to be smoothed.
            external_weights:  Weights of points in the sequence.
            bandwidth:  Smoothing bandwidth defined in terms of indices.
        
        Return value:
            Smoothed sequence.
        """
        
        y_smooth = np.empty_like(y)
        
        for i in range(len(y_smooth)):
            
            # Point indices, centred at the current point, will be used
            # as x coordinate
            x = np.arange(len(y_smooth)) - i
            
            # Compute standard weights for LOWESS
            distances = np.abs(x) / bandwidth
            weights = (1 - distances ** 3) ** 3
            np.clip(weights, 0., None, out=weights)
            
            # Include weights provided by the caller and rescale weights
            # to simplify computation of various mean values below
            weights *= external_weights
            weights /= np.sum(weights)
            
            
            # Compute smoothed value for the current point with weighted
            # least-squares fit with a linear function.  Since x
            # coordinates are centred at the current point, only need to
            # find the constant term in the linear function.
            mean_x = np.dot(weights, x)
            mean_y = np.dot(weights, y)
            mean_x2 = np.dot(weights, x ** 2)
            mean_xy = np.dot(weights, x * y)
            
            y_smooth[i] = (mean_x2 * mean_y - mean_x * mean_xy) / (mean_x2 - mean_x ** 2)
        
        return y_smooth
    
    
    def smooth(self, bandwidth=5):
        """Construct smoothed templates.
        
        Arguments:
            bandwidth:  Bandwidth for smoothing, in terms of bins.
        
        Return value:
            Tuple with smoothed templates for up and down variations.
        
        Relative up and down deviations are assumed to be symmetric in
        shape but allowed to scale differently.  The combined shape is
        smoothed using LOWESS algorithm.  The scale factors are chosen
        by minimizing chi^2 difference between smoothed and input
        templates.
        """
        
        self.smooth_averaged_deviation = self.lowess(
            0.5 * (self.up - self.down) / self.nominal, self.weights, bandwidth
        )
        
        sfUp = self._scale_factor(self.up)
        sfDown = self._scale_factor(self.down)
        
        return (
            self.nominal * (1 + sfUp * self.smooth_averaged_deviation),
            self.nominal * (1 + sfDown * self.smooth_averaged_deviation)
        )
    
    
    def _scale_factor(self, template):
        """Compute scale factor for smoothed relative deviation.
        
        Find a scale factor that gives the best match (smallest chi^2)
        between smoothed and given template.
        
        Arguments:
            template:  Template to match.
        
        Return value:
            Scale factor to be applied to smoothed relative deviation.
        """
        
        # This is result of an analytical computation
        smooth_abs_dev = self.smooth_averaged_deviation * self.nominal
        return np.sum(smooth_abs_dev * (template - self.nominal) * self.weights) / \
            np.sum(smooth_abs_dev ** 2 * self.weights)
        


if __name__ == '__main__':
    
    mpl.rc('xtick', top=True, direction='in')
    mpl.rc('ytick', right=True, direction='in')
    mpl.rc('axes', labelsize='large')
    mpl.rc('axes.formatter', limits=[-2, 4], use_mathtext=True)
    
    fig_dir = 'figSmooth'
    
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    
    
    input_file = ROOT.TFile('ttbar.root')
    output_file = ROOT.TFile('ttbar_smooth.root', 'recreate')
    templates_to_save = []
    
    binning, nominal, errors = hist_to_np(input_file.Get('TT'), full=True)
    
    
    # Copy histograms that do not need smoothing
    for key in input_file.GetListOfKeys():
        
        hist = key.ReadObj()
        name = hist.GetName()
        
        if 'MttScale' in name or 'FSR' in name or 'MassT' in name:
            continue
        
        hist.SetDirectory(output_file)
        templates_to_save.append(hist)
    
    
    for syst_name, label in [
        ('MttScale', 'Exp. scale in $m_{t\\bar t}$'),
        ('FSR', '$\\alpha_s$ in FSR'), ('MassT', '$m_t$')
    ]:
        # Smooth variations
        up = hist_to_np(input_file.Get('TT_{}Up'.format(syst_name)))
        down = hist_to_np(input_file.Get('TT_{}Down'.format(syst_name)))
        
        smoother = Smoother(nominal, up, down, 1 / errors ** 2)
        up_smooth, down_smooth = smoother.smooth()
        
        
        # Create histograms with smoothed templates
        for direction, smooth_template in [
            ('Up', up_smooth), ('Down', down_smooth)
        ]:
            hist_smooth = ROOT.TH1D(
                'TT_{}{}'.format(syst_name, direction), '',
                len(binning) - 1, binning
            )
            
            for bin in range(1, hist_smooth.GetNbinsX() + 1):
                hist_smooth.SetBinContent(bin, smooth_template[bin - 1])
            
            hist_smooth.SetDirectory(output_file)
            templates_to_save.append(hist_smooth)
        
        
        # Convert to relative deviations and plot them
        up = up / nominal - 1
        down = down / nominal - 1
        up_smooth = up_smooth / nominal - 1
        down_smooth = down_smooth / nominal - 1
        
        fig = plt.figure()
        axes = fig.add_subplot(111)
        
        axes.hist(
            binning[:-1], bins=binning, weights=up * 100, histtype='step',
            color='#a8d2f0', label='Up, input'
        )
        axes.hist(
            binning[:-1], bins=binning, weights=up_smooth * 100, histtype='step',
            color='#1f77b4', label='Up, smoothed'
        )
        axes.hist(
            binning[:-1], bins=binning, weights=down * 100, histtype='step',
            color='#ffc999', label='Down, input'
        )
        axes.hist(
            binning[:-1], bins=binning, weights=down_smooth * 100, histtype='step',
            color='#ff7f0e', label='Down, smoothed'
        )
        
        axes.axhline(0., color='black', lw=0.8, ls='dashed')
        
        axes.margins(x=0.)
        axes.set_ylim(-7.5, 7.5)
        
        axes.legend()
        axes.set_xlabel('$m_{t\\bar t}$ [GeV]')
        axes.set_ylabel('Relative deviation from nominal [%]')
        
        axes.text(
            0.5, 1.005, label, size='large',
            ha='center', va='bottom', transform=axes.transAxes
        )
        
        fig.savefig(os.path.join(fig_dir, syst_name + '.pdf'))
        plt.close(fig)
    
    
    output_file.Write()
    output_file.Close()
    input_file.Close()
