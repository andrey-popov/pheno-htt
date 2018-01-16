#include <DelphesReader.hpp>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>


DelphesReader::DelphesReader():
    bfEvent(nullptr),
    bfElectrons(nullptr), bfMuons(nullptr), bfJets(nullptr), bfMETs(nullptr),
    jetPtThreshold(20.), jetEtaThreshold(2.4)
{}


DelphesReader::~DelphesReader()
{
    delete bfEvent;
    delete bfElectrons;
    delete bfMuons;
    delete bfJets;
    delete bfMETs;
}


void DelphesReader::BeginFile(TFile *inputFile)
{
    // Set up reading of Delphes tree. It is important that all buffers are initialized with null
    //pointers (which instructs TTree to allocate memory for them).
    tree = dynamic_cast<TTree *>(inputFile->Get("Delphes"));
    
    numEvents = tree->GetEntries();
    iEvent = 0;
    
    tree->SetBranchStatus("*", false);
    
    for (auto const &mask: {"Event.Weight", "Electron.*", "Muon.*", "Jet.*", "MissingET.*"})
        tree->SetBranchStatus(mask, true);
    
    
    tree->SetBranchAddress("Event", &bfEvent);
    tree->SetBranchAddress("Electron", &bfElectrons);
    tree->SetBranchAddress("Muon", &bfMuons);
    tree->SetBranchAddress("Jet", &bfJets);
    tree->SetBranchAddress("MissingET", &bfMETs);
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


double DelphesReader::GetWeight() const
{
    return dynamic_cast<Event *>(bfEvent->At(0))->Weight;
}


Plugin::EventOutcome DelphesReader::ProcessEventToOutcome()
{
    // Check if there are events remaining in the input tree
    if (iEvent == numEvents)
        return Plugin::EventOutcome::NoEvents;
    
    
    // Read the next event
    electrons.clear();
    muons.clear();
    jets.clear();
    
    tree->GetEntry(iEvent);
    ++iEvent;
    
    
    // Copy objects from collections into vectors to avoid dealing with TCloneArrays. In case of
    //jets only save those that pass the kinematic selection.
    for (int i = 0; i < bfElectrons->GetEntries(); ++i)
        electrons.emplace_back(*dynamic_cast<Electron *>(bfElectrons->At(i)));
    
    for (int i = 0; i < bfMuons->GetEntries(); ++i)
        muons.emplace_back(*dynamic_cast<Muon *>(bfMuons->At(i)));
    
    for (int i = 0; i < bfJets->GetEntries(); ++i)
    {
        Jet *j = dynamic_cast<Jet *>(bfJets->At(i));
        
        if (j.PT > jetPtThreshold and std::abs(j.Eta) < jetEtaThreshold)
            jets.emplace_back(*j);
    }
    
    
    // Make sure collections are ordered in pt
    auto comp = [](auto const &c1, auto const &c2){return (c1.PT > c2.PT);};
    
    std::sort(electrons.begin(), electrons.end(), comp);
    std::sort(muons.begin(), muons.end(), comp);
    std::sort(jets.begin(), jets.end(), comp);
    
    
    return Plugin::EventOutcome::Ok;
}
