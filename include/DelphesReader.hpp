#pragma once

#include <DelphesReaderBase.hpp>


class TClonesArray;
class TTree;



/**
 * \class DelphesReader
 * 
 * A plugin to read Delphes files
 */
class DelphesReader: public DelphesReaderBase
{
public:
    using DelphesReaderBase::ReadOptions;
    
public:
    /**
     * Constructor from a bit mask
     * 
     * The bit mask specifies the data to read. Standard reconstructed objects and nominal event
     * weight are always read, and the mask allows to request additional data. Flags are defined by
     * enumeration ReadOptions.
     */
    DelphesReader(unsigned readOptions = 0);
    
    virtual ~DelphesReader();
    
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
    
    /**
     * Returns particles from the LHE record
     * 
     * Collection is not empty only if reading these data has been requested explicitly.
     */
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
    
    TClonesArray *bfElectrons;
    std::vector<Electron> electrons;
    
    TClonesArray *bfMuons;
    std::vector<Muon> muons;
    
    TClonesArray *bfJets;
    std::vector<Jet> jets;
    
    TClonesArray *bfMETs;
    
    TClonesArray *bfLHEParticles;
    std::vector<GenParticle> lheParticles;
};
