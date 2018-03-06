#include <DelphesReader.hpp>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>


DelphesReader::DelphesReader(unsigned readOptions):
    DelphesReaderBase(),
    bfEvents(nullptr),
    bfElectrons(nullptr), bfMuons(nullptr), bfJets(nullptr), bfMETs(nullptr),
    bfLHEParticles(nullptr)
{
    readLHEParticles = ((readOptions & LHE_PARTICLES) == LHE_PARTICLES);
}


DelphesReader::~DelphesReader()
{
    delete bfEvents;
    delete bfElectrons;
    delete bfMuons;
    delete bfJets;
    delete bfMETs;
    
    delete bfLHEParticles;
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
    
    if (readLHEParticles)
        tree->SetBranchStatus("ParticleLHEF.*", true);
    
    
    tree->SetBranchAddress("Event", &bfEvents);
    tree->SetBranchAddress("Electron", &bfElectrons);
    tree->SetBranchAddress("Muon", &bfMuons);
    tree->SetBranchAddress("Jet", &bfJets);
    tree->SetBranchAddress("MissingET", &bfMETs);
    
    if (readLHEParticles)
        tree->SetBranchAddress("ParticleLHEF", &bfLHEParticles);
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


std::vector<GenParticle> const &DelphesReader::GetLHEParticles() const
{
    if (not readLHEParticles)
        throw std::runtime_error("DelphesReader::GetLHEParticles: Reading of LHE particles "
          "has not been requested");
    
    return lheParticles;
}


MissingET const &DelphesReader::GetMissPt() const
{
    return *dynamic_cast<MissingET *>(bfMETs->At(0));
}


double DelphesReader::GetWeight() const
{
    return dynamic_cast<HepMCEvent *>(bfEvents->At(0))->Weight;
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
    
    lheParticles.clear();
    
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
    
    if (readLHEParticles)
    {
        for (int i = 0; i < bfLHEParticles->GetEntries(); ++i)
            lheParticles.emplace_back(*dynamic_cast<GenParticle *>(bfLHEParticles->At(i)));
    }
    
    
    // Make sure collections are ordered in pt
    auto comp = [](auto const &c1, auto const &c2){return (c1.PT > c2.PT);};
    
    std::sort(electrons.begin(), electrons.end(), comp);
    std::sort(muons.begin(), muons.end(), comp);
    std::sort(jets.begin(), jets.end(), comp);
    
    
    return Plugin::EventOutcome::Ok;
}
