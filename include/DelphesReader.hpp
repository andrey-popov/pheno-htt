#pragma once

#include <Plugin.hpp>

#include <classes/DelphesClasses.h>

#include <vector>


class TClonesArray;
class TTree;



/**
 * \class DelphesReader
 * 
 * A plugin to read Delphes files
 */
class DelphesReader: public Plugin
{
public:
    /// Flags to request reading of additional data
    enum ReadOptions
    {
        LHE_PARTICLES = 0x1
    };
    
public:
    /**
     * Constructor from a bit mask
     * 
     * The bit mask specifies the data to read. Standard reconstructed objects and nominal event
     * weight are always read, and the mask allows to request additional data. Flags are defined by
     * enumeration ReadOptions.
     */
    DelphesReader(unsigned readOptions = 0);
    
    ~DelphesReader();
    
public:
    /// Sets up reading of Delphes tree
    virtual void BeginFile(TFile *inputFile) override;
    
    /// Returns collection of electrons
    std::vector<Electron> const &GetElectrons() const;
    
    /// Returns collection of muons
    std::vector<Muon> const &GetMuons() const;
    
    /**
     * Returns collection of jets
     * 
     * Only jets that meet a kinematic selection are included in the collection.
     */
    std::vector<Jet> const &GetJets() const;
    
    /**
     * Returns particles from the LHE record
     * 
     * Collection is not empty only if reading these data has been requested explicitly.
     */
    std::vector<GenParticle> const &GetLHEParticles() const;
    
    /// Returns missing pt
    MissingET const &GetMissPt() const;
    
    /// Returns nominal per-event weight
    double GetWeight() const;
    
    /// Reads next event from the input file
    virtual EventOutcome ProcessEventToOutcome() override;
    
private:
    /// Non-owning pointer to Delphes tree
    TTree *tree;
    
    /// Total number of events in the tree and index of the current event
    unsigned long long numEvents, iEvent;
    
    /// Flag showing whether LHE particles should be read
    bool readLHEParticles;
    
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
    
    /// Kinematic selection applied to jets
    double jetPtThreshold, jetEtaThreshold;
};
