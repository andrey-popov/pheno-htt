#include <TTRecoPerf.hpp>

#include <TVector2.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>


TTRecoPerf::TTRecoPerf(DelphesReader const *reader_, TTReco const *ttReco_,
  std::string const &outFileName):
    reader(reader_), ttReco(ttReco_),
    outputFile(outFileName.c_str(), "recreate"),
    profBias("Bias", ";m_{tt}^{true} [GeV];Relative bias in reconstructed m_{tt}",
      12, 350., 1000.),
    profEfficiency("Efficiency", ";m_{tt}^{true} [GeV];Eff. of identification of all jets",
      12, 350., 1000.),
    nVisited(0), nTargeted(0), nReconstructable(0)
{
    profBias.SetDirectory(&outputFile);
    profEfficiency.SetDirectory(&outputFile);
}


TTRecoPerf::~TTRecoPerf()
{
    outputFile.Write();
    outputFile.Close();
}


void TTRecoPerf::PrintCounts() const
{
    std::cout << "Event counts in TTRecoPerf\n";
    std::cout << "Visited:                " << nVisited << '\n';
    std::cout << "Targeted decays at LHE: " << nTargeted << '\n';
    std::cout << "Reconstructable events: " << nReconstructable << '\n';
}


Jet const *TTRecoPerf::MatchJet(GenParticle const *p, double maxDR) const
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


bool TTRecoPerf::ProcessEvent()
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
                throw std::runtime_error("TTRecoPerf::ProcessEvent: Found more than two "
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
    
    ++nTargeted;
    
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
    Jet const *jetBLep = MatchJet(bLep, 0.4);
    Jet const *jetBHad = MatchJet(bHad, 0.4);
    Jet const *jetQ1 = MatchJet(q1, 0.4);
    Jet const *jetQ2 = MatchJet(q2, 0.4);
    
    if (jetBLep == nullptr or jetBHad == nullptr or jetQ1 == nullptr or jetQ2 == nullptr)
        return false;
    
    if (jetBLep == jetBHad or jetBLep == jetQ1 or jetBLep == jetQ2 or jetBHad == jetQ1 or
      jetBHad == jetQ2 or jetQ1 == jetQ2)
        return false;
    
    
    // Also require that jets matched to the b quarks are b-tagged
    if (jetBLep->BTag != 1 or jetBHad->BTag != 1)
        return false;
    
    ++nReconstructable;
    
    
    // Reorder jets matched to the light-flavour quarks since the ordering might be different from
    //the one at the quark level
    if (jetQ1->PT < jetQ2->PT)
        std::swap(jetQ1, jetQ2);
    
    
    // Evaluate performance
    double const trueMtt = (particles.at(bLep->M1).P4() + particles.at(bHad->M1).P4()).M();
    double const recoMtt = (ttReco->GetTopLepP4() + ttReco->GetTopHadP4()).M();
    double const weight = reader->GetWeight();
    
    profBias.Fill(trueMtt, recoMtt / trueMtt - 1., weight);
    
    bool const matched = (jetBLep == &ttReco->GetJet(TTReco::DecayJet::bTopLep) and
      jetBHad == &ttReco->GetJet(TTReco::DecayJet::bTopHad) and
      jetQ1 == &ttReco->GetJet(TTReco::DecayJet::q1TopHad) and
      jetQ2 == &ttReco->GetJet(TTReco::DecayJet::q2TopHad));
    
    profEfficiency.Fill(trueMtt, double(matched), weight);
    
    
    return true;
}
