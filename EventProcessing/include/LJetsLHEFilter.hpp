#pragma once

#include <AnalysisPlugin.hpp>

class DelphesReaderBase;


/**
 * \class LJetsLHEFilter
 * 
 * Selects events with l+jets LHE final state. Tau leptons are not allowed.
 */
class LJetsLHEFilter: public AnalysisPlugin
{
public:
    LJetsLHEFilter(DelphesReaderBase const *reader);
    
private:
    /// Applies event selection
    virtual bool ProcessEvent() override;
    
private:
    /// Non-owning pointer to reader plugin
    DelphesReaderBase const *reader;
};
