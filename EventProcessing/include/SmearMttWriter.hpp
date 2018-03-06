#pragma once

#include <AnalysisPlugin.hpp>

#include <DelphesReaderBase.hpp>

#include <TRandom3.h>


class TTree;


/**
 * \class SmearMttWriter
 * 
 * A plugin that stores in a tree the parton-level mass of the tt system and a smeared value that
 * mimics reconstruction effects.
 */
class SmearMttWriter: public AnalysisPlugin
{
public:
    SmearMttWriter(DelphesReaderBase const *reader, double resolution);
    
public:
    virtual void BeginFile(TFile *) override;
    
private:
    virtual bool ProcessEvent() override;
    
private:
    /// Non-owning pointer to reader plugin
    DelphesReaderBase const *reader;
    
    /// Relative resolution to smear mtt
    double resolution;
    
    /// Random-number generator used for smearing
    TRandom3 rGen;
    
    /// Non-owning pointer to output tree
    TTree *outTree;
    
    // Buffers to fill the tree
    Float_t bfWeight;
    Float_t bfPartonMassTT, bfMassTT;
};
