#!/usr/bin/env python

"""Computes significance and CLs for hMSSM.

Performs a scan over the (mA, tan(beta)) plane.
"""

import argparse
import itertools
import json

import numpy as np
from scipy.integrate import trapz

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT.RooStats import HistFactory

import lhapdf

import hmssm
from signalmtt import SignalMtt


def build_signal_templates(mA, tanBeta, binning, resolution):
    """Construct templates for signal.
    
    Templates are split into contributions with positive and negative
    bin contents, and the sign of the latter ones is flipped.  As a
    result, they should be scaled with negative signal strength in the
    statistical model.  Construct also templates representing a
    variation of the renormalization scale by factor 2.  Produced
    histograms have trivial equidistant binning.
    
    Arguments:
        mA, tanBeta:  Parameters of hMSSM.
        binning:  Binning in reconstructed mtt.
        resolution:  Relative resolution in reconstructed mtt.
    
    Return value:
        A dictionary with signal templates.
    """
    
    templates = {}
    
    partonCalc = hmssm.PartonXSec(mA, tanBeta, 'hMSSM_13TeV.root')
    signalMtt = SignalMtt(partonCalc, resolution=resolution)
    
    for muRScaleFactor, postfix in [
        (1., ''), (2., '_RenormScaleUp'), (0.5, '_RenormScaleDown')
    ]:
        signalMtt.muRScaleFactor = muRScaleFactor
        signalMtt.build_xsec_grid()
        
        
        # In each bin given by the binning, integrate the differential
        # cross section
        xSecIntegral = np.empty(len(binning) - 1)
        nSamples = 10
        
        for iBin in range(len(xSecIntegral)):
            
            x = np.linspace(binning[iBin], binning[iBin + 1], num=nSamples)
            y = np.empty_like(x)
            
            for iSample in range(nSamples):
                y[iSample] = signalMtt.xsec(x[iSample])
            
            xSecIntegral[iBin] = trapz(y, x)
        
        
        # Put results into histograms, separately for positive and
        # negative counts
        histPos = ROOT.TH1D(
            'SgnPos' + postfix, '',
            len(binning) - 1, 0., len(binning) - 1
        )
        histNeg = histPos.Clone('SgnNeg' + postfix)
        
        for iBin in range(len(xSecIntegral)):
            value = xSecIntegral[iBin]
            
            if value >= 0.:
                histPos.SetBinContent(iBin + 1, value)
            else:
                histNeg.SetBinContent(iBin + 1, value)
            
            histPos.SetBinError(iBin + 1, 0.)
            histNeg.SetBinError(iBin + 1, 0.)
        
        histNeg.Scale(-1.)
        
        for hist in [histPos, histNeg]:
            hist.SetDirectory(None)
            templates[hist.GetName()] = hist
    
    return templates


def build_model(signalTemplates, bkgFile, lumi):
    """Build RooFit workspace with statistical model.
    
    Arguments:
        signalTemplates:  Dictionary with signal templates.
        bkgFile:  ROOT file with templates for SM tt background.
        lumi:  Target luminosity, 1/fb.
    
    Return value:
        RooFit workspace with statistical model.
    """
    
    measurement = HistFactory.Measurement('Measurement')
    
    measurement.AddPOI('r')
    measurement.AddPreprocessFunction('negR', '-r', 'r[1,0,10]')
    measurement.AddConstantParam('Lumi')
    measurement.SetLumi(lumi)
    
    
    channel = HistFactory.Channel('Channel')
    
    sgnPos = HistFactory.Sample('SgnPos')
    sgnPos.SetHisto(signalTemplates['SgnPos'])
    sgnPos.AddNormFactor('r', 1., 0., 10., True)
    sgnPos.SetNormalizeByTheory(True)
    
    syst = HistFactory.HistoSys('RenormScale')
    syst.SetHistoHigh(signalTemplates['SgnPos_RenormScaleUp'])
    syst.SetHistoLow(signalTemplates['SgnPos_RenormScaleDown'])
    sgnPos.AddHistoSys(syst)
    
    sgnNeg = HistFactory.Sample('SgnNeg')
    sgnNeg.SetHisto(signalTemplates['SgnNeg'])
    sgnNeg.AddNormFactor('negR', -1., -10., 0., True)
    sgnNeg.SetNormalizeByTheory(True)
    
    syst = HistFactory.HistoSys('RenormScale')
    syst.SetHistoHigh(signalTemplates['SgnNeg_RenormScaleUp'])
    syst.SetHistoLow(signalTemplates['SgnNeg_RenormScaleDown'])
    sgnNeg.AddHistoSys(syst)
    
    bkg = HistFactory.Sample('TT')
    bkg.SetHisto(bkgFile.Get('TT'))
    bkg.AddOverallSys('TTRate', 0.9, 1.1)
    
    systNames = ['MttScale', 'RenormScale', 'FactorScale', 'FSR', 'MassT', 'PDFAlphaS']
    systNames += ['PDF{}'.format(i) for i in range(1, 31)]
    
    for systName in systNames:
        syst = HistFactory.HistoSys(systName)
        syst.SetHistoHigh(bkgFile.Get('TT_{}Up'.format(systName)))
        syst.SetHistoLow(bkgFile.Get('TT_{}Down'.format(systName)))
        bkg.AddHistoSys(syst)
    
    
    for sample in [sgnPos, sgnNeg, bkg]:
        channel.AddSample(sample)
    
    measurement.AddChannel(channel)
    
    
    # An s+b Asimov data set will be created automatically.  Add also a
    # b-only data set.
    bOnlyData = HistFactory.Asimov('asimovData_bOnly')
    bOnlyData.SetParamValue('r', 0.)
    measurement.AddAsimovDataset(bOnlyData)
    
    workspace = HistFactory.HistoToWorkspaceFactoryFast.MakeCombinedModel(measurement)
    return workspace


def compute_significance(workspace):
    """Compute significance on s+b Asimov data set.
    
    Quantifies the incompatibility between the s+b Asimov data set and
    post-fit SM expectation.
    
    Arguments:
        workspace:  RooFit workspace defining the statistical model.
    
    Return value:
        Expected significance.
    """
    
    fullModel = workspace.obj('ModelConfig')
    
    # Create a snapshot that fixes the POI in the full model
    poi = fullModel.GetParametersOfInterest().first()
    poi.setVal(1.)
    fullModel.SetSnapshot(ROOT.RooArgSet(poi))
    
    # Create the b-only model by fixing POI to zero in the full model
    bOnlyModel = fullModel.Clone(fullModel.GetName() + '_bOnly')
    poi.setVal(0.)
    bOnlyModel.SetSnapshot(ROOT.RooArgSet(poi))
    
    
    # Asymptotic calculator for the significance with Asimov s+b data
    # set as the data
    calc = ROOT.RooStats.AsymptoticCalculator(
        workspace.data('asimovData'), fullModel, bOnlyModel
    )
    
    # Configure the calculator to use the q0 statistics from [1], which
    # is intended for discovery
    # [1] https://arxiv.org/abs/1007.1727
    calc.SetOneSidedDiscovery(True)
    
    
    # Perform the test.  The Asimov data set provided to the constructor
    # of calc is fitted.  Because of the data set used, the "observed"
    # significance coincides with median expected one.
    testResult = calc.GetHypoTest()
    
    return testResult.Significance()


def compute_CLs(workspace):
    """Compute CLs on b-only Asimov data set.
    
    Arguments:
        workspace:  RooFit workspace defining the statistical model.
    
    Return value:
        CLs value.
    """
    
    fullModel = workspace.obj('ModelConfig')
    
    # Create a snapshot that fixes the POI in the full model
    poi = fullModel.GetParametersOfInterest().first()
    poi.setVal(1.)
    fullModel.SetSnapshot(ROOT.RooArgSet(poi))
    
    # Create the b-only model by fixing POI to zero in the full model
    bOnlyModel = fullModel.Clone(fullModel.GetName() + '_bOnly')
    poi.setVal(0.)
    bOnlyModel.SetSnapshot(ROOT.RooArgSet(poi))
    
    
    # Asymptotic calculator for the significance with Asimov b-only data
    # set as the data.  For an exclusion the null hypothesis is that
    # r >= 1 and the alternative hypothesis is r = 0.
    calc = ROOT.RooStats.AsymptoticCalculator(
        workspace.data('asimovData_bOnly'), bOnlyModel, fullModel
    )
    
    # Configure the calculator to use the q_mu tilde statistics from
    # [1], which is intended for upper limits with additional condition
    # r > 0.
    # [1] https://arxiv.org/abs/1007.1727
    calc.SetOneSided(True)
    calc.SetQTilde(True)
    
    
    # Perform the test.  The Asimov data set provided to the constructor
    # of calc is fitted.  Because of the data set used, the "observed"
    # CLs coincides with median expected one.
    testResult = calc.GetHypoTest()
    testResult.SetBackgroundAsAlt(True)
    
    return testResult.CLs()


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
        '-o', '--output', default='scan.csv',
        help='Name for output CSV file'
    )
    args = argParser.parse_args()
    
    ROOT.Math.MinimizerOptions.SetDefaultMinimizer('Minuit2', 'Minimize')
    ROOT.RooStats.UseNLLOffset(True)
    lhapdf.setVerbosity(0)
    
    
    with open(args.binning) as f:
        binning = np.array(json.load(f), dtype=np.float64)
    
    bkgFile = ROOT.TFile(args.bkg)
    
    
    results = []
    
    for mA, tanBeta in itertools.product(
        np.arange(350, 1001, 25), np.arange(0.75, 5.1, 0.25)
    ):
        signalTemplates = build_signal_templates(mA, tanBeta, binning, args.resolution)
        workspace = build_model(signalTemplates, bkgFile, args.lumi * 1e3)
        significance = compute_significance(workspace)
        cls = compute_CLs(workspace)
        
        print('\033[1;34mResults for mA = {:g}, tan(beta) = {:g}:'.format(mA, tanBeta))
        print('  Significance: {}\n  CLs:  {}\033[0m'.format(significance, cls))
        results.append((mA, tanBeta, significance, cls))
    
    
    bkgFile.Close()
    
    outputFile = open(args.output, 'w')
    outputFile.write('#mA,tanBeta,Significance,CLs\n')
    
    for entry in results:
        outputFile.write('{},{},{},{}\n'.format(*entry))
    
    outputFile.close()
