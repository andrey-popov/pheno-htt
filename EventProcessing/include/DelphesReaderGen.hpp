#pragma once

#include <DelphesReaderBase.hpp>


class TClonesArray;
class TTree;



/**
 * \class DelphesReaderGen
 * 
 * A plugin to read generator-level information from Delphes files. It is meant to be a drop-in
 * replacement of class DelphesReader when the reconstruction has not been performed.
 */
class DelphesReaderGen: public DelphesReaderBase
{
public:
    DelphesReaderGen();
    virtual ~DelphesReaderGen();
    
public:
    /// Sets up reading of Delphes tree
    virtual void BeginFile(TFile *inputFile) override;
    
    /// Returns collection of electrons
    virtual std::vector<Electron> const &GetElectrons() const override;
    
    /// Returns collection of muons
    virtual std::vector<Muon> const &GetMuons() const override;
    
    /**
     * Returns collection of jets
     * 
     * Only jets that meet a kinematic selection are included in the collection.
     */
    virtual std::vector<Jet> const &GetJets() const override;
    
    /// Returns particles from the LHE record
    virtual std::vector<GenParticle> const &GetLHEParticles() const override;
    
    /// Returns missing pt
    virtual MissingET const &GetMissPt() const override;
    
    /// Returns nominal per-event weight
    virtual double GetWeight() const override;
    
    /// Reads next event from the input file
    virtual EventOutcome ProcessEventToOutcome() override;
    
private:
    /// Non-owning pointer to Delphes tree
    TTree *tree;
    
    /// Total number of events in the tree and index of the current event
    unsigned long long numEvents, iEvent;
    
    /// Buffer to read global generator-level information about an event
    TClonesArray *bfEvents;
    
    TClonesArray *bfLHEParticles;
    std::vector<GenParticle> lheParticles;
    
    TClonesArray *bfJets;
    std::vector<Jet> jets;
    
    TClonesArray *bfMETs;
    
    std::vector<Electron> electrons;
    std::vector<Muon> muons;
};
