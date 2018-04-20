#include <DelphesReader.hpp>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>


DelphesReader::DelphesReader():
    DelphesReaderBase(),
    bfElectrons(nullptr), bfMuons(nullptr), bfJets(nullptr), bfMETs(nullptr)
{}


DelphesReader::~DelphesReader()
{
    delete bfElectrons;
    delete bfMuons;
    delete bfJets;
    delete bfMETs;
}


std::vector<Electron> const &DelphesReader::GetElectrons() const
{
    return electrons;
}


std::vector<Muon> const &DelphesReader::GetMuons() const
{
    return muons;
}


std::vector<Jet> const &DelphesReader::GetJets() const
{
    return jets;
}


MissingET const &DelphesReader::GetMissPt() const
{
    return *dynamic_cast<MissingET *>(bfMETs->At(0));
}


void DelphesReader::ReadEvent()
{
    // Read objects provided by the base class
    DelphesReaderBase::ReadEvent();
    
    
    // Copy objects from collections into vectors to avoid dealing with TCloneArrays. In case of
    //jets only save those that pass the kinematic selection.
    electrons.clear();
    muons.clear();
    jets.clear();
    
    for (int i = 0; i < bfElectrons->GetEntries(); ++i)
        electrons.emplace_back(*dynamic_cast<Electron *>(bfElectrons->At(i)));
    
    for (int i = 0; i < bfMuons->GetEntries(); ++i)
        muons.emplace_back(*dynamic_cast<Muon *>(bfMuons->At(i)));
    
    for (int i = 0; i < bfJets->GetEntries(); ++i)
    {
        Jet *j = dynamic_cast<Jet *>(bfJets->At(i));
        
        if (j->PT < jetPtThreshold or std::abs(j->Eta) > jetEtaThreshold)
            continue;
        
        // There is a memory leak in TRefArray [1], and class Jet contains members of this type.
        //The leak would occur when the collection of jets is sorted below. To prevent it, clear
        //the arrays completely.
        //[1] https://sft.its.cern.ch/jira/browse/ROOT-7589
        j->Constituents.Delete();
        j->Particles.Delete();
        
        jets.emplace_back(std::move(*j));
    }
    
    
    // Make sure collections are ordered in pt
    auto comp = [](auto const &c1, auto const &c2){return (c1.PT > c2.PT);};
    
    std::sort(electrons.begin(), electrons.end(), comp);
    std::sort(muons.begin(), muons.end(), comp);
    std::sort(jets.begin(), jets.end(), comp);
}


void DelphesReader::SetupBuffers()
{
    // Setup buffers for the base class
    DelphesReaderBase::SetupBuffers();
    
    
    for (auto const &mask: {"Electron.*", "Muon.*", "Jet.*", "MissingET.*"})
        tree->SetBranchStatus(mask, true);    
    
    tree->SetBranchAddress("Electron", &bfElectrons);
    tree->SetBranchAddress("Muon", &bfMuons);
    tree->SetBranchAddress("Jet", &bfJets);
    tree->SetBranchAddress("MissingET", &bfMETs);
}
