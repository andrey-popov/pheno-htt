#include <LJetsSelection.hpp>

#include <cmath>


LJetsSelection::LJetsSelection(DelphesReader *reader_):
    reader(reader_),
    ptEleTight(30.), ptEleLoose(10.),
    ptMuTight(30.), ptMuLoose(10.),
    mtWThreshold(0.)
{}


TLorentzVector const &LJetsSelection::GetLeptonP4() const
{
    return p4TightLepton;
}


double LJetsSelection::GetMtW() const
{
    return mtW;
}


bool LJetsSelection::ProcessEvent()
{
    // Count tight and loose leptons and save four-momentum of the only tight lepton
    unsigned nTight = 0, nLoose = 0;
    
    for (auto const &e: reader->GetElectrons())
    {
        if (e.PT < ptEleLoose or std::abs(e.Eta) > 2.5)
            continue;
        
        ++nLoose;
        
        if (e.PT > ptEleTight)
        {
            ++nTight;
            p4TightLepton = e.P4();
        }
    }
    
    for (auto const &mu: reader->GetMuons())
    {
        if (mu.PT < ptMuLoose or std::abs(mu.Eta) > 2.4)
            continue;
        
        ++nLoose;
        
        if (mu.PT > ptMuTight)
        {
            ++nTight;
            p4TightLepton = mu.P4();
        }
    }
    
    if (nTight != 1 or nLoose != 1)
        return false;
    
    
    // Count jets. Selection on them has already been applied by the reader.
    auto const &jets = reader->GetJets();
    
    if (jets.size() < 4)
        return false;
    
    unsigned nTags = 0;
    
    for (auto const &j: jets)
    {
        if (j.BTag == 1)
            ++nTags;
    }
    
    if (nTags < 2)
        return false;
    
    
    // Compute m_T(W) and apply selection on it
    auto const &met = reader->GetMissPt();
    mtW = std::sqrt(2 * p4TightLepton.Pt() * met.MET *
      (1 - std::cos(p4TightLepton.Phi() - met.Phi)));
    
    if (mtW < mtWThreshold)
        return false;
    
    
    // If the workflow reaches this point, the event is accepted
    return true;
}
