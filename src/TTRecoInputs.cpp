#include <TTRecoInputs.hpp>

#include <LJetsSelection.hpp>

#include <TVector2.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>


TTRecoInputs::TTRecoInputs(DelphesReader *reader_, LJetsSelection *selector_,
  std::string const &outFileName):
    reader(reader_), selector(selector_),
    outputFile(outFileName.c_str(), "recreate"),
    histNeutrinoDist("NeutrinoDist", ";|#vec{p}_{T}^{miss} - #vec{p}_{T}^{#nu}| [GeV];Events",
      100, 0., 200.),
    histMassesHad("MassesHad", ";m_{t} [GeV];m_{W} [GeV];Events", 200, 50., 250., 150, 0., 150.),
    nuReco(173., 80.419002),  // Masses from param_card
    nVisited(0), nTargetLHE(0), nReconstructable(0), nFilled(0)
{
    histNeutrinoDist.SetDirectory(&outputFile);
    histMassesHad.SetDirectory(&outputFile);
}


TTRecoInputs::~TTRecoInputs()
{
    outputFile.Write();
    outputFile.Close();
}


void TTRecoInputs::PrintCounts() const
{
    std::cout << "Event counts in TTRecoInputs\n";
    std::cout << "Visited:                " << nVisited << '\n';
    std::cout << "Targeted decays at LHE: " << nTargetLHE << '\n';
    std::cout << "Reconstructable events: " << nReconstructable << '\n';
    std::cout << "Filled in histograms:   " << nFilled << std::endl;
}


Jet const *TTRecoInputs::MatchJet(GenParticle const *p, double maxDR) const
{
    auto const &jets = reader->GetJets();
    Jet const *match = nullptr;
    double minDR2 = std::pow(maxDR, 2);
    
    for (auto const &j: jets)
    {
        double const dR2 = std::pow(p->Eta - j.Eta, 2) +
          std::pow(TVector2::Phi_mpi_pi(p->Phi - j.Phi), 2);
        
        if (dR2 < minDR2)
        {
            match = &j;
            minDR2 = dR2;
        }
    }
    
    return match;
}


bool TTRecoInputs::ProcessEvent()
{
    ++nVisited;
    auto const &particles = reader->GetLHEParticles();
    
    // Select events with targeted decays at the LHE level and identify b quarks and light-flavour
    //quarks from decays of W bosons
    unsigned nLep = 0, nTau = 0, nB = 0, nQ = 0;
    std::array<GenParticle const *, 2> bQuarks, lightQuarks;
    
    for (auto const &p: particles)
    {
        int const absPdgID = std::abs(p.PID);
        
        if (absPdgID == 11 or absPdgID == 13)
            ++nLep;
        else if (absPdgID == 15)
            ++nTau;
        else if (absPdgID == 5 and std::abs(particles.at(p.M1).PID) == 6)
        {
            if (nB == 2)
                throw std::runtime_error("TTRecoInputs::ProcessEvent: Found more than two "
                  "b quarks.");
            
            bQuarks[nB] = &p;
            ++nB;
        }
        else if (absPdgID <= 4 and absPdgID != 0 and p.M1 != -1 and
          std::abs(particles.at(p.M1).PID) == 24)
        {
            if (nQ == 2)
            {
                // This cannot be the targeted decay
                return false;
            }
            
            lightQuarks[nQ] = &p;
            ++nQ;
        }
    }
    
    if (nLep != 1 or nTau > 0)
        return false;
    
    ++nTargetLHE;
    
    assert(nB == 2);
    assert(nQ == 2);
    
    
    // Order light-flavour quarks by pt and distinguish b quarks from semileptonic and hadronic
    //decays
    GenParticle const *q1 = lightQuarks[0], *q2 = lightQuarks[1];
    
    if (q1->PT < q2->PT)
        std::swap(q1, q2);
    
    GenParticle const *bLep = bQuarks[0], *bHad = bQuarks[1];
    
    if (bLep->M1 == particles.at(q1->M1).M1)
        std::swap(bLep, bHad);
    
    assert(bHad->M1 == particles.at(q1->M1).M1);
    
    
    // Check if the quarks can be matched to reconstructed jets
    Jet const *jetBLep = MatchJet(bLep);
    Jet const *jetBHad = MatchJet(bHad);
    Jet const *jetQ1 = MatchJet(q1);
    Jet const *jetQ2 = MatchJet(q2);
    
    if (jetBLep == nullptr or jetBHad == nullptr or jetQ1 == nullptr or jetQ2 == nullptr)
        return false;
    
    if (jetBLep == jetBHad or jetBLep == jetQ1 or jetBLep == jetQ2 or jetBHad == jetQ1 or
      jetBHad == jetQ2 or jetQ1 == jetQ2)
        return false;
    
    
    // Also require that jets matched to the b quarks are b-tagged
    if (jetBLep->BTag != 1 or jetBHad->BTag != 1)
        return false;
    
    ++nReconstructable;
    
    
    // Perform neutrino reconstruction
    nuReco.Reconstruct(selector->GetLeptonP4(), jetBLep->P4(), reader->GetMissPt().P4());
    
    if (nuReco.RecoStatus() != 0)
        return false;
    
    
    // Finally, fill the histograms
    double const weight = reader->GetWeight();
    histNeutrinoDist.Fill(nuReco.GetCompatibility(), weight);
    
    TLorentzVector const p4WHad = jetQ1->P4() + jetQ2->P4();
    histMassesHad.Fill((p4WHad + jetBHad->P4()).M(), p4WHad.M(), weight);
    
    ++nFilled;
    
    
    return true;
}
