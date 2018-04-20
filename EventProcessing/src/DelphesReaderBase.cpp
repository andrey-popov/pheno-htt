#include <DelphesReaderBase.hpp>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>


DelphesReaderBase::DelphesReaderBase(double jetPtThreshold_, double jetEtaThreshold_):
    jetPtThreshold(jetPtThreshold_), jetEtaThreshold(jetEtaThreshold_),
    bfEvents(nullptr), bfLHEParticles(nullptr)
{}


DelphesReaderBase::~DelphesReaderBase()
{
    delete bfEvents;
    delete bfLHEParticles;
}


void DelphesReaderBase::BeginFile(TFile *inputFile)
{
    // Set up reading of Delphes tree. It is important that all buffers are initialized with null
    //pointers (which instructs TTree to allocate memory for them).
    tree = dynamic_cast<TTree *>(inputFile->Get("Delphes"));
    
    numEvents = tree->GetEntries();
    iEvent = 0;
    
    tree->SetBranchStatus("*", false);
    SetupBuffers();
}


std::vector<GenParticle> const &DelphesReaderBase::GetLHEParticles() const
{
    return lheParticles;
}


double DelphesReaderBase::GetWeight() const
{
    return dynamic_cast<HepMCEvent *>(bfEvents->At(0))->Weight;
}


Plugin::EventOutcome DelphesReaderBase::ProcessEventToOutcome()
{
    // Check if there are events remaining in the input tree
    if (iEvent == numEvents)
        return Plugin::EventOutcome::NoEvents;
    
    
    // Read the next event
    tree->GetEntry(iEvent);
    ++iEvent;
    
    ReadEvent();
    
    return Plugin::EventOutcome::Ok;
}


void DelphesReaderBase::ReadEvent()
{
    // Read LHE particles. Copy objects from collections into vectors to avoid dealing with
    //TCloneArrays.
    lheParticles.clear();
    
    for (int i = 0; i < bfLHEParticles->GetEntries(); ++i)
        lheParticles.emplace_back(*dynamic_cast<GenParticle *>(bfLHEParticles->At(i)));
}


void DelphesReaderBase::SetupBuffers()
{
    for (auto const &mask: {"Event.Weight", "ParticleLHEF.*"})
        tree->SetBranchStatus(mask, true);
    
    tree->SetBranchAddress("Event", &bfEvents);
    tree->SetBranchAddress("ParticleLHEF", &bfLHEParticles);
}
