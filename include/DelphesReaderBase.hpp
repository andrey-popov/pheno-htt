#pragma once

#include <Plugin.hpp>

#include <classes/DelphesClasses.h>

#include <vector>


/**
 * \class DelphesReaderBase
 * 
 * An abstract base class for a plugin that reads Delphes files
 */
class DelphesReaderBase: public Plugin
{
public:
    /// Constructor that defines kinematic selection applied to jets
    DelphesReaderBase(double jetPtThreshold = 20., double jetEtaThreshold = 2.4);
    
    virtual ~DelphesReaderBase() = default;
    
public:
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
    
    /**
     * Returns particles from the LHE record
     * 
     * Collection is not empty only if reading these data has been requested explicitly.
     */
    virtual std::vector<GenParticle> const &GetLHEParticles() const = 0;
    
    /// Returns missing pt
    virtual MissingET const &GetMissPt() const = 0;
    
    /// Returns nominal per-event weight
    virtual double GetWeight() const = 0;
    
protected:
    /// Kinematic selection applied to jets
    double jetPtThreshold, jetEtaThreshold;
};
