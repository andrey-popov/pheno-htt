#include <VarWriter.hpp>

#include <Processor.hpp>


VarWriter::VarWriter(DelphesReaderBase const *reader_, TTReco const *ttReco_):
    reader(reader_), ttReco(ttReco_),
    storePartonLevel(false)
{}


void VarWriter::BeginFile(TFile *)
{
    outTree = processor->Book<TTree>("", "Vars", "Observables computed for tt system");
    
    outTree->Branch("Weight", &bfWeight);
    outTree->Branch("PtTopLep", &bfPtTopLep);
    outTree->Branch("PtTopHad", &bfPtTopHad);
    outTree->Branch("MassTT", &bfMassTT);
    
    if (storePartonLevel)
        outTree->Branch("PartonMassTT", &bfPartonMassTT);
}


void VarWriter::StorePartonLevel(bool on)
{
    storePartonLevel = on;
}


bool VarWriter::ProcessEvent()
{
    bfWeight = reader->GetWeight();
    
    TLorentzVector const p4TopLep = ttReco->GetTopLepP4();
    TLorentzVector const p4TopHad = ttReco->GetTopHadP4();
    
    bfPtTopLep = p4TopLep.Pt();
    bfPtTopHad = p4TopHad.Pt();
    bfMassTT = (p4TopLep + p4TopHad).M();
    
    
    if (storePartonLevel)
    {
        auto const &particles = reader->GetLHEParticles();
        
        TLorentzVector p4TT;
        
        for (auto const &p: particles)
        {
            if (std::abs(p.PID) == 6)
                p4TT += p.P4();
        }
        
        bfPartonMassTT = p4TT.M();
    }
    
    
    outTree->Fill();
    return true;
}
