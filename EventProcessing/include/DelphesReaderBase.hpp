#pragma once

#include <Plugin.hpp>

#include <classes/DelphesClasses.h>

#include <vector>

class TClonesArray;
class TTree;


/**
 * \class DelphesReaderBase
 * 
 * An abstract base class for a plugin that reads Delphes files
 * 
 * It reads certain generator-level information and provides an interface to access reconstructed
 * objects. Derived classes are provided with hooks SetupBuffers and ReadEvent to allow reading
 * additional collections from a Delphes file.
 */
class DelphesReaderBase: public Plugin
{
public:
    /// Constructor that defines kinematic selection applied to jets
    DelphesReaderBase(double jetPtThreshold = 20., double jetEtaThreshold = 2.4);
    
    virtual ~DelphesReaderBase();
    
public:
    /// Sets up reading of Delphes tree
    virtual void BeginFile(TFile *inputFile) override final;
    
    /// Returns collection of electrons
    virtual std::vector<Electron> const &GetElectrons() const = 0;
    
    /// Returns collection of muons
    virtual std::vector<Muon> const &GetMuons() const = 0;
    
    /**
     * Returns collection of jets
     * 
     * Only jets that meet a kinematic selection are included in the collection.
     */
    virtual std::vector<Jet> const &GetJets() const = 0;
    
    /// Returns particles from the LHE record
    std::vector<GenParticle> const &GetLHEParticles() const;
    
    /**
     * Returns event weights from the LHE record
     * 
     * Only available if reading of LHE weights has been requested with SetReadLHEWeights
     * beforehand.
     */
    std::vector<LHEFWeight> const &GetLHEWeights() const;
    
    /// Returns missing pt
    virtual MissingET const &GetMissPt() const = 0;
    
    /// Returns nominal per-event weight
    double GetWeight() const;
    
    /// Reads next event from the input file
    virtual EventOutcome ProcessEventToOutcome() override final;
    
    /// Requests reading of LHE weights
    void SetReadLHEWeights(bool on = true);
    
protected:
    /**
     * A hook to read content of the current event
     * 
     * A derived class must use this hook to read additional collections from the event. It should
     * also call the method of the base class to read collections provided by DelphesReaderBase.
     */
    virtual void ReadEvent();
    
    /**
     * A hood to set up buffers to read branches of the Delphes trees
     * 
     * A derived class must use this hook to enable branches it requires (by default all branches
     * are disabled) and assign buffers for them. It should also call the method of the base class
     * to allow reading of collections provided by DelphesReaderBase.
     */
    virtual void SetupBuffers();
    
protected:
    /// Kinematic selection applied to jets
    double jetPtThreshold, jetEtaThreshold;
    
    /// Non-owning pointer to Delphes tree
    TTree *tree;
    
    /// Total number of events in the tree and index of the current event
    unsigned long long numEvents, iEvent;
    
    /// Buffer to read global generator-level information about an event
    TClonesArray *bfEvents;
    
    /// Buffer to read LHE particles
    TClonesArray *bfLHEParticles;
    std::vector<GenParticle> lheParticles;
    
    /// Buffer to read LHE weights
    TClonesArray *bfLHEWeights;
    std::vector<LHEFWeight> lheWeights;
    
    /// Indicates whether LHE weights should be read
    bool readLHEWeights;
};
