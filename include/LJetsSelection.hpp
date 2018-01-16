#pragma once

#include <AnalysisPlugin.hpp>

#include <DelphesReader.hpp>

#include <TLorentzVector.h>


/**
 * \class LJetsSelection
 * 
 * A plugin to implement tt -> l+jets selection
 * 
 * Event is required to contain exactly one tight electron or muon, no additional leptons, at least
 * four jets, and at least two of them must be b-tagged. Additionally, the value of m_T(W) must be
 * above a threshold.
 */
class LJetsSelection: public AnalysisPlugin
{
public:
    /// Constructor from a pointer to reader plugin
    LJetsSelection(DelphesReader *reader);
    
public:
    /// Returns four-momentum of the only tight lepton in the current accepted event
    TLorentzVector const &GetLeptonP4() const;
    
    /// Returns value of m_T(W) computed for the current accepted event
    double GetMtW() const;
    
    /// Applies event selection
    virtual bool ProcessEvent() override;
    
private:
    /// Non-owning pointer to reader plugin
    DelphesReader *reader;
    
    /// Thresholds for electrons
    double ptEleTight, ptEleLoose;
    
    /// Thresholds for jets
    double ptMuTight, ptMuLoose;
    
    /// Four-momentum of the only tight lepton allowed in the event
    TLorentzVector p4TightLepton;
    
    /// Value of the m_T(W) variable
    double mtW;
    
    /// Minimal threshold for m_T(W) for the event selection
    double mtWThreshold;
};
