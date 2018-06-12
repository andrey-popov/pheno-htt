#!/usr/bin/env python

"""Removes binning information from histograms.

Reads all histograms in input ROOT file and copies their content into
new histograms with equidistant binning with unit bin width.  The new
histograms are then saved in output ROOT file.  This effectively removes
binning information from the histograms.  This procedure is needed for
HistFactor as it seems to misbehave when histograms have a
non-equidistant binning [1].

[1] https://root-forum.cern.ch/t/binned-ml-fit-with-histfactory/28867
"""

import argparse

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


if __name__ == '__main__':
    
    argParser = argparse.ArgumentParser(epilog=__doc__)
    argParser.add_argument('inputFile', help='Input ROOT file')
    argParser.add_argument('outputFile', help='Output ROOT file')
    args = argParser.parse_args()
    
    
    inputFile = ROOT.TFile(args.inputFile)
    outputFile = ROOT.TFile(args.outputFile, 'recreate')
    histsToStore = []
    
    for key in inputFile.GetListOfKeys():
        
        # Only 1D histograms are considered
        if not key.GetClassName().startswith('TH1'):
            continue
        
        sourceHist = key.ReadObj()
        
        nBins = sourceHist.GetNbinsX()
        transformedHist = ROOT.TH1D(sourceHist.GetName(), '', nBins, 0., nBins)
        transformedHist.SetDirectory(outputFile)
        
        # Copy content of the source histogram.  Under- and overflow
        # bins are ignored.
        for bin in range(1, nBins + 1):
            transformedHist.SetBinContent(bin, sourceHist.GetBinContent(bin))
            transformedHist.SetBinError(bin, sourceHist.GetBinError(bin))
        
        
        # Prevent garbage collection from removing this histogram
        histsToStore.append(transformedHist)
    
    outputFile.Write()
    outputFile.Close()
    inputFile.Close()
