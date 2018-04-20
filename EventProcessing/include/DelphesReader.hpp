#pragma once

#include <DelphesReaderBase.hpp>


class TClonesArray;



/**
 * \class DelphesReader
 * 
 * A plugin to read Delphes files
 */
class DelphesReader: public DelphesReaderBase
{
public:
    DelphesReader();
    
    virtual ~DelphesReader();
    
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
    TClonesArray *bfElectrons;
    std::vector<Electron> electrons;
    
    TClonesArray *bfMuons;
    std::vector<Muon> muons;
    
    TClonesArray *bfJets;
    std::vector<Jet> jets;
    
    TClonesArray *bfMETs;
};
