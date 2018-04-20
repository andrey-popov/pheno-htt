#include <DelphesReaderGen.hpp>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector.h>

#include <boost/container/small_vector.hpp>

#include <algorithm>
#include <cmath>
#include <functional>


template<typename T, typename C>
bool CheckOverlap(T const &obj, C const &collection, double dR = 0.4)
{
    for (auto const &element: collection)
    {
        // Make a reference_wrapper regardless of whether elements of the collection are bare types
        //or instances of reference_wrapper so that the two cases can be treated uniformly
        auto const elementRef = std::cref(element);
        
        double const dR2 = std::pow(obj.Eta - elementRef.get().Eta, 2) +
          std::pow(TVector2::Phi_mpi_pi(obj.Phi - elementRef.get().Phi), 2);
        
        if (dR2 < dR * dR)
            return true;
    }
    
    return false;
}


DelphesReaderGen::DelphesReaderGen():
    DelphesReaderBase(),
    bfJets(nullptr), bfMETs(nullptr)
{}


DelphesReaderGen::~DelphesReaderGen()
{
    delete bfJets;
    delete bfMETs;
}


std::vector<Electron> const &DelphesReaderGen::GetElectrons() const
{
    return electrons;
}


std::vector<Muon> const &DelphesReaderGen::GetMuons() const
{
    return muons;
}


std::vector<Jet> const &DelphesReaderGen::GetJets() const
{
    return jets;
}


MissingET const &DelphesReaderGen::GetMissPt() const
{
    return *dynamic_cast<MissingET *>(bfMETs->At(0));
}


void DelphesReaderGen::ReadEvent()
{
    // Read objects provided by the base class
    DelphesReaderBase::ReadEvent();
    
    
    // Fill collections of leptons with muons and electrons from LHE. Only some fields are set.
    electrons.clear();
    muons.clear();
    
    for (auto const &p: lheParticles)
    {
        if (std::abs(p.PID) == 11)
        {
            Electron e;
            
            e.PT = p.PT;
            e.Eta = p.Eta;
            e.Phi = p.Phi;
            e.Charge = (p.PID > 0) ? -1 : +1;
            
            e.T = e.EhadOverEem = e.IsolationVar = e.IsolationVarRhoCorr = e.SumPtCharged = \
              e.SumPtNeutral = e.SumPtChargedPU = e.SumPt = 0.;
            
            electrons.emplace_back(std::move(e));
        }
        else if (std::abs(p.PID) == 13)
        {
            Muon mu;
            
            mu.PT = p.PT;
            mu.Eta = p.Eta;
            mu.Phi = p.Phi;
            mu.Charge = (p.PID > 0) ? -1 : +1;
            
            mu.T = mu.IsolationVar = mu.IsolationVarRhoCorr = mu.SumPtCharged = mu.SumPtNeutral = \
              mu.SumPtChargedPU = mu.SumPt = 0.;
            
            muons.emplace_back(std::move(mu));
        }
    }
    
    
    // Copy jets from the read buffer. Only keep those that pass the kinematic selection and do not
    //overlap with the leptons.
    jets.clear();
    
    for (int i = 0; i < bfJets->GetEntries(); ++i)
    {
        Jet *j = dynamic_cast<Jet *>(bfJets->At(i));
        
        if (j->PT < jetPtThreshold or std::abs(j->Eta) > jetEtaThreshold)
            continue;
        
        if (CheckOverlap(*j, electrons) or CheckOverlap(*j, muons))
            continue;
        
        // There is a memory leak in TRefArray [1], and class Jet contains members of this type.
        //The leak would occur when the collection of jets is sorted below. To prevent it, clear
        //the arrays completely.
        //[1] https://sft.its.cern.ch/jira/browse/ROOT-7589
        j->Constituents.Delete();
        j->Particles.Delete();
        
        jets.emplace_back(std::move(*j));
    }
    
    
    // Emulate b-tagging by matching jets to b quarks
    boost::container::small_vector<std::reference_wrapper<GenParticle const>, 2> bQuarks;
    
    for (auto const &p: lheParticles)
    {
        if (std::abs(p.PID) == 5)
            bQuarks.emplace_back(p);
    }
    
    for (auto &j: jets)
    {
        if (CheckOverlap(j, bQuarks))
            j.BTag = 1;
    }
    
    
    // Make sure collections are ordered in pt
    auto comp = [](auto const &c1, auto const &c2){return (c1.PT > c2.PT);};
    
    std::sort(electrons.begin(), electrons.end(), comp);
    std::sort(muons.begin(), muons.end(), comp);
    std::sort(jets.begin(), jets.end(), comp);
}


void DelphesReaderGen::SetupBuffers()
{
    // Setup buffers for the base class
    DelphesReaderBase::SetupBuffers();
    
    
    for (auto const &mask: {"GenJet.*", "GenMissingET.*"})
        tree->SetBranchStatus(mask, true);
    
    tree->SetBranchAddress("GenJet", &bfJets);
    tree->SetBranchAddress("GenMissingET", &bfMETs);
}
