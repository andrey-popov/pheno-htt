#include <SmearMttWriter.hpp>

#include <Processor.hpp>

#include <TTree.h>


SmearMttWriter::SmearMttWriter(DelphesReaderBase const *reader_, double resolution_):
    reader(reader_), resolution(resolution_),
    rGen(0)
{}


void SmearMttWriter::BeginFile(TFile *)
{
    outTree = processor->Book<TTree>("", "Vars", "Mass of tt system");
    
    outTree->Branch("Weight", &bfWeight);
    outTree->Branch("PartonMassTT", &bfPartonMassTT);
    outTree->Branch("MassTT", &bfMassTT);
}


bool SmearMttWriter::ProcessEvent()
{
    bfWeight = reader->GetWeight();
    
    
    // Compute parton-level mass
    auto const &particles = reader->GetLHEParticles();
    TLorentzVector p4TT;
    
    for (auto const &p: particles)
    {
        if (std::abs(p.PID) == 6)
            p4TT += p.P4();
    }
    
    bfPartonMassTT = p4TT.M();
    
    
    // Smear the mass
    bfMassTT = rGen.Gaus(bfPartonMassTT, bfPartonMassTT * resolution);
    
    
    outTree->Fill();
    return true;
}
