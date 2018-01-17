#pragma once

#include <AnalysisPlugin.hpp>

#include <DelphesReader.hpp>
#include <TTReco.hpp>

#include <string>

#include <TFile.h>
#include <TProfile.h>


/**
 * \class TTRecoPerf
 * 
 * An auxiliary plugin to evaluate performance of tt reconstruction
 * 
 * This plugin selects reconstructable tt -> l+jets events and fills two TProfiles, both in bins of
 * parton-level mtt. The first one is filled with relative bias in reconstructed mtt (and bin
 * uncertainties quantify the resolution). The second is filled with 1 or 0 depending on whether
 * the four jets have been identified correctly; the mean values then give the reconstruction
 * efficiency. The profiles are computed with all input files and stored in a signle ROOT file.
 */
class TTRecoPerf: public AnalysisPlugin
{
public:
    /// Constructor with a name for the output file
    TTRecoPerf(DelphesReader const *reader, TTReco const *ttReco, std::string const &outFileName);
    
    ~TTRecoPerf();
    
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
    DelphesReader const *reader;
    
    /// Non-owning pointer to a plugin that performs tt reconstruction
    TTReco const *ttReco;
    
    /// Output ROOT file
    TFile outputFile;
    
    /// Profiles of mtt bias and efficiency of jet identification
    TProfile profBias, profEfficiency;
    
    /// Event counters for sanity checks
    unsigned long long nVisited, nTargeted, nReconstructable;
};
