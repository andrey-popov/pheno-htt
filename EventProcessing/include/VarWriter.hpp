#pragma once

#include <AnalysisPlugin.hpp>

#include <DelphesReaderBase.hpp>
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
    VarWriter(DelphesReaderBase const *reader, TTReco const *ttReco);
    
public:
    virtual void BeginFile(TFile *) override;
    
    /// Requests storing of some parton-level observables
    void StorePartonLevel(bool on = true);
    
private:
    virtual bool ProcessEvent() override;
    
private:
    /// Non-owning pointer to reader plugin
    DelphesReaderBase const *reader;
    
    /// Non-owning pointer to plugin that performs tt reconstruction
    TTReco const *ttReco;
    
    /// Flag showing whether parton-level variables should be stored
    bool storePartonLevel;
    
    /// Non-owning pointer to output tree
    TTree *outTree;
    
    // Buffers to fill the tree
    Float_t bfWeight;
    Float_t bfPtTopLep, bfPtTopHad, bfMassTT;
    Float_t bfPartonMassTT;
};
