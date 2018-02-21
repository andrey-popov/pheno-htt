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
    /// Flags to request reading of additional data
    enum ReadOptions
    {
        LHE_PARTICLES = 0x1
    };
    
public:
    /**
     * Constructor from a bit mask
     * 
     * The mask allows to request reading of additional data (to be implemented in a derived
     * class). Supported flags are defined by enumeration ReadOptions.
     */
    DelphesReaderBase(unsigned readOptions = 0);
    
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
    /// Flag showing whether LHE particles should be read
    bool readLHEParticles;
    
    /// Kinematic selection applied to jets
    double jetPtThreshold, jetEtaThreshold;
};
