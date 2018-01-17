#pragma once

#include <AnalysisPlugin.hpp>

#include <DelphesReader.hpp>
#include <TTReco.hpp>

#include <TTree.h>


/**
 * \class VarWriter
 * 
 * A pluging to compute various observables and store them in ROOT trees
 */
class VarWriter: public AnalysisPlugin
{
public:
    VarWriter(DelphesReader const *reader, TTReco const *ttReco);
    
public:
    virtual void BeginFile(TFile *) override;
    
private:
    virtual bool ProcessEvent() override;
    
private:
    /// Non-owning pointer to reader plugin
    DelphesReader const *reader;
    
    /// Non-owning pointer to plugin that performs tt reconstruction
    TTReco const *ttReco;
    
    /// Non-owning pointer to output tree
    TTree *outTree;
    
    // Buffers to fill the tree
    Float_t bfWeight;
    Float_t bfPtTopLep, bfPtTopHad, bfMassTT;
};
