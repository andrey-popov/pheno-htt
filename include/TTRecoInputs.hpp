#pragma once

#include <AnalysisPlugin.hpp>

#include <DelphesReader.hpp>
#include <NuReco.hpp>

#include <string>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>


class DelphesReader;
class LJetsSelection;


/**
 * \class TTRecoInputs
 * 
 * An auxiliary plugin to produce inputs needed for reconstruction of tt -> l+jets decays
 * 
 * This plugin selects tt -> l+jets events that can in principle be fully reconstructed and fills
 * histograms of Euclidian distance between measured missing pt and pt of reconstructed neutrino,
 * for the semileptonic leg of the decay, and for the hadronic leg it fills a 2D histogram of
 * masses of the top quark and the W boson computed from reconstructed jets. The histograms are
 * combined over all input files and stored in a ROOT file.
 */
class TTRecoInputs: public AnalysisPlugin
{
public:
    /// Constructor with a name for the output file
    TTRecoInputs(DelphesReader *reader, LJetsSelection *selector, std::string const &outFileName);
    
    ~TTRecoInputs();
    
public:
    /// Prints a summary of event counts
    void PrintCounts() const;
    
private:
    /**
     * Matches a reconstructed jet to the given GenParticle
     * 
     * Matching is done by dR. If not match with dR less than the given threshold is found, a null
     * pointer is returned.
     */
    Jet const *MatchJet(GenParticle const *p, double maxDR = 0.2) const;
    
    /// Selects suitable events and fills the histograms
    virtual bool ProcessEvent() override;
    
private:
    /// Non-owning pointer to reader plugin
    DelphesReader *reader;
    
    /// Non-owning pointer to a plugin that performs event selection
    LJetsSelection *selector;
    
    /// Output ROOT file
    TFile outputFile;
    
    /// Histogram of Euclidian distance between missing pt and pt of reconstructed neutrino
    TH1D histNeutrinoDist;
    
    /// Histogram of (mt, mW)
    TH2D histMassesHad;
    
    /// An object to perform neutrino reconstruction
    NuReco nuReco;
    
    /// Event counters for sanity checks
    unsigned long long nVisited, nTargetLHE, nReconstructable, nFilled;
};
