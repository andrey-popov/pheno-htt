#pragma once

#include <AnalysisPlugin.hpp>

#include <DelphesReaderBase.hpp>

#include <TRandom3.h>

#include <string>
#include <vector>


class TH1;


/**
 * \class SystMttHists
 * 
 * A plugin to create histograms of smeared parton-level mtt
 * 
 * This plugin computes parton-level mass of the tt system and applies to it a Gaussian smearing
 * with a relative resolution specified by user. It stores in the output file a nominal histogram
 * and histograms for several systematic variations, all with given binning. For the systematic
 * variations, it creates a histogram for each provided alternative LHE weight and also a pair of
 * histograms for a scale variation in the smeared mtt (which approximates a variation due to a
 * jet momentum scale uncertainty).
 */
class SystMttHists: public AnalysisPlugin
{
public:
    /**
     * Constructor from a reader, binning for mtt histograms, relative resolution, and scale
     * variation
    */
    SystMttHists(DelphesReaderBase const *reader, std::vector<double> const &binning,
      double resolution, double scaleVariation = 0.01);
    
public:
    virtual void BeginFile(TFile *) override;
    
private:
    /// Fills nominal histogram and histograms for systematic variations
    virtual bool ProcessEvent() override;
    
private:
    /// Non-owning pointer to reader plugin
    DelphesReaderBase const *reader;
    
    /// Random-number generator used for smearing
    TRandom3 rGen;
    
    /// Binning for mtt histograms
    std::vector<double> binning;
    
    /// Relative resolution for mtt
    double resolution;
    
    /// Scale variation for mtt
    double scaleVariation;
    
    /// Non-owning pointer to nominal histogram
    TH1 *histNominal;
    
    /// Non-owning pointers to histograms with varied mtt scale
    TH1 *histScaleUp, *histScaleDown;
    
    /// Non-owning pointers to histograms with alternative LHE weights
    std::vector<TH1 *> histAltWeights;
    
    /// Indicates whether hitograms with alternative weights have been booked
    bool histAltWeightsBooked;
};
