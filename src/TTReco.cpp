#include <TTReco.hpp>

#include <LJetsSelection.hpp>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>
#include <sstream>
#include <stdexcept>


TTReco::TTReco(DelphesReader *reader_, LJetsSelection *selector_,
  std::string const &likelihoodFileName):
    reader(reader_), selector(selector_),
    nuReco(173., 80.419002),  // Masses from param_card
    minPt(0.), maxAbsEta(std::numeric_limits<double>::infinity())
{
    // Read histograms that define likelihood function
    TFile likelihoodFile(likelihoodFileName.c_str());
    
    if (likelihoodFile.IsZombie())
    {
        std::ostringstream message;
        message << "TTReco::TTReco: Failed to open file \"" << likelihoodFileName <<
          "\" for reading.";
        throw std::runtime_error(message.str());
    }
    
    likelihoodNuDist = dynamic_cast<TH1 *>(likelihoodFile.Get("NeutrinoDist"));
    likelihoodMassHad = dynamic_cast<TH2 *>(likelihoodFile.Get("MassesHad"));
    
    likelihoodNuDist->SetDirectory(nullptr);
    likelihoodMassHad->SetDirectory(nullptr);
    
    likelihoodFile.Close();
    
    
    // Normalize the histograms
    likelihoodNuDist->Scale(1. / likelihoodNuDist->Integral(), "width");
    likelihoodMassHad->Scale(1. / likelihoodMassHad->Integral(), "width");
}


Jet const &TTReco::GetJet(DecayJet type) const
{
    Jet const *jet = nullptr;
    
    switch (type)
    {
        case DecayJet::bTopLep:
            jet = bTopLep;
            break;
        
        case DecayJet::bTopHad:
            jet = bTopHad;
            break;
        
        case DecayJet::q1TopHad:
            jet = q1TopHad;
            break;
        
        case DecayJet::q2TopHad:
            jet = q2TopHad;
            break;
        
        default:
            // This must never happen
            throw std::runtime_error("TTReco::GetJet: Unhandled jet type.");
    }
    
    
    if (not jet)
        throw std::runtime_error("TTReco::GetJet: Requested jet is not available. "
          "This probably means that reconstruction for the current event has been aborted.");
    
    return *jet;
}


TLorentzVector const &TTReco::GetLeptonP4() const
{
    return selector->GetLeptonP4();
}


TLorentzVector const &TTReco::GetNeutrinoP4() const
{
    return nuReco.GetSolution();
}


double TTReco::GetRank() const
{
    return highestRank;
}


unsigned TTReco::GetRecoStatus() const
{
    return recoStatus;
}


TLorentzVector TTReco::GetTopLepP4() const
{
    return GetLeptonP4() + GetNeutrinoP4() + GetJet(DecayJet::bTopLep).P4();
}


TLorentzVector TTReco::GetTopHadP4() const
{
    return GetJet(DecayJet::bTopHad).P4() + GetJet(DecayJet::q1TopHad).P4() +
      GetJet(DecayJet::q2TopHad).P4();
}


void TTReco::SetJetSelection(double minPt_, double maxAbsEta_)
{
    minPt = minPt_;
    maxAbsEta = maxAbsEta_;
}


bool TTReco::ProcessEvent()
{
    // Reset data describing the current-best interpretation
    highestRank = -std::numeric_limits<double>::infinity();
    bTopLep = bTopHad = q1TopHad = q2TopHad = nullptr;
    
    
    // Apply kinematic selection to jets
    auto const &jets = reader->GetJets();
    selectedJetIndices.clear();
    
    for (unsigned i = 0; i < jets.size(); ++i)
    {
        if (std::abs(jets[i].Eta) > maxAbsEta)
            continue;
        
        if (jets[i].PT < minPt)
            break;  // The jet collection is ordered in pt
        
        selectedJetIndices.push_back(i);
    }
    
    unsigned const nSelectedJets = selectedJetIndices.size();
    
    
    // Do not attempt reconstruction if there is not enough jets
    if (nSelectedJets < 4)
    {
        recoStatus = 1;
        return false;
    }
    
    
    // Loop over all possible ways of jet assignment to find the best one
    for (unsigned iiBTopLepCand = 0; iiBTopLepCand < nSelectedJets; ++iiBTopLepCand)
    {
        // Jets matched to b quarks must be b-tagged
        if (jets.at(selectedJetIndices.at(iiBTopLepCand)).BTag != 1)
            continue;
        
        
        // Reconstruct neutrino
        nuReco.Reconstruct(GetLeptonP4(), jets.at(selectedJetIndices.at(iiBTopLepCand)).P4(),
          reader->GetMissPt().P4());
        
        if (nuReco.RecoStatus() != 0)
            continue;
        
        
        // Compute the likelihood for neutrino reconstruction
        int bin = likelihoodNuDist->FindFixBin(nuReco.GetCompatibility());
        
        if (likelihoodNuDist->IsBinOverflow(bin))
            continue;
        
        double const curLLNu = std::log(likelihoodNuDist->GetBinContent(bin));
        
        
        // Check permutations for the hadronic leg of the decay
        for (unsigned iiBTopHadCand = 0; iiBTopHadCand < nSelectedJets; ++iiBTopHadCand)
        {
            if (iiBTopLepCand == iiBTopHadCand)
                continue;
            
            // Jets matched to b quarks must be b-tagged
            if (jets.at(selectedJetIndices.at(iiBTopHadCand)).BTag != 1)
                continue;
            
            
            for (unsigned iiQ1TopHadCand = 0; iiQ1TopHadCand < nSelectedJets; ++iiQ1TopHadCand)
            {
                if (iiQ1TopHadCand == iiBTopLepCand or iiQ1TopHadCand == iiBTopHadCand)
                    continue;
                
                // When looping for the subleading light-flavour jet, take into account that the
                //collection is still ordered in jet pt
                for (unsigned iiQ2TopHadCand = iiQ1TopHadCand + 1; iiQ2TopHadCand < nSelectedJets;
                  ++iiQ2TopHadCand)
                {
                    if (iiQ2TopHadCand == iiBTopLepCand or iiQ2TopHadCand == iiBTopHadCand)
                        continue;
                    
                    // An interpretation has been constructed. Compute the full likelihood for it.
                    TLorentzVector const p4W =
                      jets.at(selectedJetIndices.at(iiQ1TopHadCand)).P4() +
                      jets.at(selectedJetIndices.at(iiQ2TopHadCand)).P4();
                    int bin = likelihoodMassHad->FindFixBin(
                      (p4W + jets.at(selectedJetIndices.at(iiBTopHadCand)).P4()).M(), p4W.M());
                    
                    if (likelihoodMassHad->IsBinOverflow(bin))
                        continue;
                    
                    double const rank = curLLNu + std::log(likelihoodMassHad->GetBinContent(bin));
                    
                    if (rank > highestRank)
                    {
                        highestRank = rank;
                        
                        bTopLep = &jets.at(selectedJetIndices.at(iiBTopLepCand));
                        bTopHad = &jets.at(selectedJetIndices.at(iiBTopHadCand));
                        q1TopHad = &jets.at(selectedJetIndices.at(iiQ1TopHadCand));
                        q2TopHad = &jets.at(selectedJetIndices.at(iiQ2TopHadCand));
                    }
                }
            }
        }
    }
    
    if (highestRank == -std::numeric_limits<double>::infinity())
    {
        recoStatus = 2;
        return false;
    }
    
    
    recoStatus = 0;
    return true;
}
