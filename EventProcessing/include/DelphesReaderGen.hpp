#pragma once

#include <DelphesReaderBase.hpp>


class TClonesArray;



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
    
    /// Returns missing pt
    virtual MissingET const &GetMissPt() const override;
    
protected:
    /// Reads additional collections from the current event
    virtual void ReadEvent() override;
    
    /// Sets up buffers to read branches of Delphes tree with additional collections
    virtual void SetupBuffers() override;
    
private:
    TClonesArray *bfJets;
    std::vector<Jet> jets;
    
    TClonesArray *bfMETs;
    
    std::vector<Electron> electrons;
    std::vector<Muon> muons;
};
