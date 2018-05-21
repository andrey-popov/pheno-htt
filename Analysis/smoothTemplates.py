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
    
    nBins = hist.GetNbinsX()
    content = np.empty(nBins)
    
    for bin in range(1, nBins + 1):
        content[bin - 1] = hist.GetBinContent(bin)
    
    if not full:
        return content
    else:
        binning = np.empty(nBins + 1)
        errors = np.empty(nBins)
        
        for bin in range(1, nBins + 1):
            binning[bin - 1] = hist.GetBinLowEdge(bin)
            errors[bin - 1] = hist.GetBinError(bin)
        
        binning[-1] = hist.GetBinLowEdge(nBins + 1)
        
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
        
        self.smoothAveragedDeviation = None
    
    
    @staticmethod
    def lowess(y, externalWeights, bandwidth):
        """Peform LOWESS smoothing.
        
        Arguments:
            y:  Input sequence to be smoothed.
            externalWeights:  Weights of points in the sequence.
            bandwidth:  Smoothing bandwidth defined in terms of indices.
        
        Return value:
            Smoothed sequence.
        """
        
        ySmooth = np.empty_like(y)
        
        for i in range(len(ySmooth)):
            
            # Point indices, centred at the current point, will be used
            # as x coordinate
            x = np.arange(len(ySmooth)) - i
            
            # Compute standard weights for LOWESS
            distances = np.abs(x) / bandwidth
            weights = (1 - distances ** 3) ** 3
            np.clip(weights, 0., None, out=weights)
            
            # Include weights provided by the caller and rescale weights
            # to simplify computation of various mean values below
            weights *= externalWeights
            weights /= np.sum(weights)
            
            
            # Compute smoothed value for the current point with weighted
            # least-squares fit with a linear function.  Since x
            # coordinates are centred at the current point, only need to
            # find the constant term in the linear function.
            meanX = np.dot(weights, x)
            meanY = np.dot(weights, y)
            meanX2 = np.dot(weights, x ** 2)
            meanXY = np.dot(weights, x * y)
            
            ySmooth[i] = (meanX2 * meanY - meanX * meanXY) / (meanX2 - meanX ** 2)
        
        return ySmooth
    
    
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
        
        self.smoothAveragedDeviation = self.lowess(
            0.5 * (self.up - self.down) / self.nominal, self.weights, bandwidth
        )
        
        sfUp = self._scale_factor(self.up)
        sfDown = self._scale_factor(self.down)
        
        return (
            self.nominal * (1 + sfUp * self.smoothAveragedDeviation),
            self.nominal * (1 + sfDown * self.smoothAveragedDeviation)
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
        smoothAbsDev = self.smoothAveragedDeviation * self.nominal
        return np.sum(smoothAbsDev * (template - self.nominal) * self.weights) / \
            np.sum(smoothAbsDev ** 2 * self.weights)
        


if __name__ == '__main__':
    
    mpl.rc('xtick', top=True, direction='in')
    mpl.rc('ytick', right=True, direction='in')
    mpl.rc('axes', labelsize='large')
    mpl.rc('axes.formatter', limits=[-2, 4], use_mathtext=True)
    
    figDir = 'figSmooth'
    
    if not os.path.exists(figDir):
        os.makedirs(figDir)
    
    
    inputFile = ROOT.TFile('ttbar.root')
    outputFile = ROOT.TFile('ttbar_smooth.root', 'recreate')
    templatesToSave = []
    
    binning, nominal, errors = hist_to_np(inputFile.Get('TT'), full=True)
    
    
    # Copy histograms that do not need smoothing
    for key in inputFile.GetListOfKeys():
        
        hist = key.ReadObj()
        name = hist.GetName()
        
        if 'MttScale' in name or 'FSR' in name or 'MassT' in name:
            continue
        
        hist.SetDirectory(outputFile)
        templatesToSave.append(hist)
    
    
    for systName, label in [
        ('MttScale', 'Exp. scale in $m_{t\\bar t}$'),
        ('FSR', '$\\alpha_s$ in FSR'), ('MassT', '$m_t$')
    ]:
        # Smooth variations
        up = hist_to_np(inputFile.Get('TT_{}Up'.format(systName)))
        down = hist_to_np(inputFile.Get('TT_{}Down'.format(systName)))
        
        smoother = Smoother(nominal, up, down, 1 / errors ** 2)
        upSmooth, downSmooth = smoother.smooth()
        
        
        # Create histograms with smoothed templates
        for direction, smoothTemplate in [
            ('Up', upSmooth), ('Down', downSmooth)
        ]:
            histSmooth = ROOT.TH1D(
                'TT_{}{}'.format(systName, direction), '',
                len(binning) - 1, binning
            )
            
            for bin in range(1, histSmooth.GetNbinsX() + 1):
                histSmooth.SetBinContent(bin, smoothTemplate[bin - 1])
            
            histSmooth.SetDirectory(outputFile)
            templatesToSave.append(histSmooth)
        
        
        # Convert to relative deviations and plot them
        up = up / nominal - 1
        down = down / nominal - 1
        upSmooth = upSmooth / nominal - 1
        downSmooth = downSmooth / nominal - 1
        
        fig = plt.figure()
        axes = fig.add_subplot(111)
        
        axes.hist(
            binning[:-1], bins=binning, weights=up * 100, histtype='step',
            color='#a8d2f0', label='Up, input'
        )
        axes.hist(
            binning[:-1], bins=binning, weights=upSmooth * 100, histtype='step',
            color='#1f77b4', label='Up, smoothed'
        )
        axes.hist(
            binning[:-1], bins=binning, weights=down * 100, histtype='step',
            color='#ffc999', label='Down, input'
        )
        axes.hist(
            binning[:-1], bins=binning, weights=downSmooth * 100, histtype='step',
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
        
        fig.savefig(os.path.join(figDir, systName + '.pdf'))
        plt.close(fig)
    
    
    outputFile.Write()
    outputFile.Close()
    inputFile.Close()
