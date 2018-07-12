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
    
    def __init__(self, output_name, num_events, k_factor=1.):
        """Initializer.
        
        Arguments:
            output_name:  Name for output ROOT file.
            num_events:  Number of generated events before event
                selection.
            k_factor:  Scale factor to account for higher-order
                corrections.
        """
        
        self.output_file = ROOT.TFile(output_name, 'recreate')
        self.k_factor = k_factor
        self.num_events = num_events
        self.hists = []
    
    
    def add(self, input_file, hist_name, hist_write_name=None):
        """Add a new histogram to the output file.
        
        Arguments:
            input_file:  Opened ROOT file or a path to one.
            hist_name:  Name of histogram to read from the file.
            hist_write_name:  New name for the histogram, which will be
                used in the output file.
        
        Return value:
            Added histogram
        """
        
        if not isinstance(input_file, ROOT.TDirectoryFile):
            # Interpret this argument as a path
            input_file = ROOT.TFile(input_file)
            opened_file = True
        else:
            opened_file = False
        
        hist = input_file.Get(hist_name)
        
        if not hist:
            raise RuntimeError('Failed to read histogram "{}{}".'.format(
                input_file.GetPath(), hist_name
            ))
        
        if hist_write_name:
            hist.SetName(hist_write_name)
        
        hist.Scale(self.k_factor / self.num_events)
        
        # Associate the histogram with the output file and put it into a
        # list so that it is not deleted by garbage collector
        hist.SetDirectory(self.output_file)
        self.hists.append(hist)
        
        # Close input file if it was opened here
        if opened_file:
            input_file.Close()
        
        return hist
    
    
    def add_hist(self, hist):
        """Add fully constructed histogram.
        
        Do not apply any modifications to it.
        """
        
        hist.SetDirectory(self.output_file)
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
        
        self.output_file.Write()


if __name__ == '__main__':
    
    ROOT.gROOT.SetBatch(True)
    
    aggregator = HistAggregator('ttbar.root', 5000000, k_factor=1.6)
    
    main_input_file = ROOT.TFile('hists/ttbar.root')
    
    hist_nominal = aggregator.add(main_input_file, 'Nominal', 'TT')
    
    aggregator.add(main_input_file, 'ScaleUp', 'TT_MttScaleUp')
    aggregator.add(main_input_file, 'ScaleDown', 'TT_MttScaleDown')
    
    
    # Changes in cross section are factorized out so that the resulting
    # variations only reflect the impact on the acceptance and shapes
    hist = aggregator.add(main_input_file, 'AltWeight_ID35', 'TT_RenormScaleUp')
    hist.Scale(78.5413478374 / 66.0966981708)
    
    hist = aggregator.add(main_input_file, 'AltWeight_ID6',  'TT_RenormScaleDown')
    hist.Scale(78.5413478374 / 94.9362566747)
    
    hist = aggregator.add(main_input_file, 'AltWeight_ID25', 'TT_FactorScaleUp')
    hist.Scale(78.5413478374 / 74.2105427156)
    
    hist = aggregator.add(main_input_file, 'AltWeight_ID16', 'TT_FactorScaleDown')
    hist.Scale(78.5413478374 / 83.0940774172)
    
    
    for ipdf in range(30):
        hist_pdf_up = aggregator.add(
            main_input_file, 'AltWeight_ID{}'.format(46 + ipdf), 'TT_PDF{}Up'.format(ipdf + 1)
        )
        
        # PDF uncertainties are described with symmetric eigenvectors,
        # and only one variation encoded into the weights.  It is taken
        # as the "up" variation.  The "down" one is symmetric by
        # construction.  It is built manually.
        hist_pdf_down = hist_nominal.Clone('TT_PDF{}Down'.format(ipdf + 1))
        hist_pdf_down.Scale(2)
        hist_pdf_down.Add(hist_pdf_up, -1)
        aggregator.add_hist(hist_pdf_down)
        
    
    aggregator.add(main_input_file, 'AltWeight_ID77', 'TT_PDFAlphaSUp')
    aggregator.add(main_input_file, 'AltWeight_ID76', 'TT_PDFAlphaSDown')
    
    main_input_file.Close()
    
    aggregator.add('hists/ttbar_FSR-up.root', 'Nominal', 'TT_FSRUp')
    aggregator.add('hists/ttbar_FSR-down.root', 'Nominal', 'TT_FSRDown')
    
    aggregator.add('hists/ttbar_mt-up.root', 'Nominal', 'TT_MassTUp')
    aggregator.add('hists/ttbar_mt-down.root', 'Nominal', 'TT_MassTDown')
    
    aggregator.save()
