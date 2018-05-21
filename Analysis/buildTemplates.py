#!/usr/bin/env python

"""Constructs SM tt templates for RooStats.

All parameters, including names of input and output files, total number
of events before the event selection, and k-factor, are hard-coded.

Variations in the renormalization and factorization scales in the matrix
element are rescaled to the nominal cross section.
"""

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


class HistAggregator:
    """A class to simplify aggregation of different histograms.
    
    It copies histograms from multiple sources into a single output
    file.  The histograms can be renamed.  They are rescaled for the
    number of generated events, and a k-factor is applied.
    """
    
    def __init__(self, outputName, numEvents, kFactor=1.):
        """Initializer.
        
        Arguments:
            outputName:  Name for output ROOT file.
            numEvents:  Number of generated events before event
                selection.
            kFactor:  Scale factor to account for higher-order
                corrections.
        """
        
        self.outputFile = ROOT.TFile(outputName, 'recreate')
        self.kFactor = kFactor
        self.numEvents = numEvents
        self.hists = []
    
    
    def add(self, inputFile, histName, histWriteName=None):
        """Add a new histogram to the output file.
        
        Arguments:
            inputFile:  Opened ROOT file or a path to one.
            histName:  Name of histogram to read from the file.
            histWriteName:  New name for the histogram, which will be
                used in the output file.
        
        Return value:
            Added histogram
        """
        
        if not isinstance(inputFile, ROOT.TDirectoryFile):
            # Interpret this argument as a path
            inputFile = ROOT.TFile(inputFile)
            openedFile = True
        else:
            openedFile = False
        
        hist = inputFile.Get(histName)
        
        if not hist:
            raise RuntimeError('Failed to read histogram "{}{}".'.format(
                inputFile.GetPath(), histName
            ))
        
        if histWriteName:
            hist.SetName(histWriteName)
        
        hist.Scale(self.kFactor / self.numEvents)
        
        # Associate the histogram with the output file and put it into a
        # list so that it is not deleted by garbage collector
        hist.SetDirectory(self.outputFile)
        self.hists.append(hist)
        
        # Close input file if it was opened here
        if openedFile:
            inputFile.Close()
        
        return hist
    
    
    def add_hist(self, hist):
        """Add fully constructed histogram.
        
        Do not apply any modifications to it.
        """
        
        hist.SetDirectory(self.outputFile)
        self.hists.append(hist)
    
    
    def save(self):
        """Save histograms associated with output file.
        
        Also clean under- and overflows in all histograms.
        """
        
        for hist in self.hists:
            hist.SetBinContent(0, 0.)
            hist.SetBinError(0, 0.)
            hist.SetBinContent(hist.GetNbinsX() + 1, 0.)
            hist.SetBinError(hist.GetNbinsX() + 1, 0.)
        
        self.outputFile.Write()


if __name__ == '__main__':
    
    ROOT.gROOT.SetBatch(True)
    
    aggregator = HistAggregator('ttbar.root', 5000000, kFactor=1.6)
    
    mainInputFile = ROOT.TFile('hists/ttbar.root')
    
    histNominal = aggregator.add(mainInputFile, 'Nominal', 'TT')
    
    aggregator.add(mainInputFile, 'ScaleUp', 'TT_MttScaleUp')
    aggregator.add(mainInputFile, 'ScaleDown', 'TT_MttScaleDown')
    
    
    # Changes in cross section are factorized out so that the resulting
    # variations only reflect the impact on the acceptance and shapes
    hist = aggregator.add(mainInputFile, 'AltWeight_ID35', 'TT_RenormScaleUp')
    hist.Scale(78.5413478374 / 66.0966981708)
    
    hist = aggregator.add(mainInputFile, 'AltWeight_ID6',  'TT_RenormScaleDown')
    hist.Scale(78.5413478374 / 94.9362566747)
    
    hist = aggregator.add(mainInputFile, 'AltWeight_ID25', 'TT_FactorScaleUp')
    hist.Scale(78.5413478374 / 74.2105427156)
    
    hist = aggregator.add(mainInputFile, 'AltWeight_ID16', 'TT_FactorScaleDown')
    hist.Scale(78.5413478374 / 83.0940774172)
    
    
    for iPDF in range(30):
        histPDFUp = aggregator.add(
            mainInputFile, 'AltWeight_ID{}'.format(46 + iPDF), 'TT_PDF{}Up'.format(iPDF + 1)
        )
        
        # PDF uncertainties are described with symmetric eigenvectors,
        # and only one variation encoded into the weights.  It is taken
        # as the "up" variation.  The "down" one is symmetric by
        # construction.  It is built manually.
        histPDFDown = histNominal.Clone('TT_PDF{}Down'.format(iPDF + 1))
        histPDFDown.Scale(2)
        histPDFDown.Add(histPDFUp, -1)
        aggregator.add_hist(histPDFDown)
        
    
    aggregator.add(mainInputFile, 'AltWeight_ID77', 'TT_PDFAlphaSUp')
    aggregator.add(mainInputFile, 'AltWeight_ID76', 'TT_PDFAlphaSDown')
    
    mainInputFile.Close()
    
    aggregator.add('hists/ttbar_FSR-up.root', 'Nominal', 'TT_FSRUp')
    aggregator.add('hists/ttbar_FSR-down.root', 'Nominal', 'TT_FSRDown')
    
    aggregator.add('hists/ttbar_mt-up.root', 'Nominal', 'TT_MassTUp')
    aggregator.add('hists/ttbar_mt-down.root', 'Nominal', 'TT_MassTDown')
    
    aggregator.save()
