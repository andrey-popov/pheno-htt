#include "SystMttHists.hpp"

#include <Processor.hpp>

#include <TH1D.h>


SystMttHists::SystMttHists(DelphesReaderBase const *reader_, std::vector<double> const &binning_,
  double resolution_, double scaleVariation_):
    reader(reader_), rGen(0),
    binning(binning_), resolution(resolution_), scaleVariation(scaleVariation_)
{}


void SystMttHists::BeginFile(TFile *)
{
    histNominal = processor->Book<TH1D>("", "Nominal", "", binning.size() - 1, binning.data());
    histScaleUp = processor->Book<TH1D>("", "ScaleUp", "", binning.size() - 1, binning.data());
    histScaleDown = processor->Book<TH1D>("", "ScaleDown", "", binning.size() - 1, binning.data());
    
    histAltWeightsBooked = false;
    histAltWeights.clear();
}


bool SystMttHists::ProcessEvent()
{
    // Compute parton-level mass
    auto const &particles = reader->GetLHEParticles();
    TLorentzVector p4TT;
    
    for (auto const &p: particles)
    {
        if (std::abs(p.PID) == 6)
            p4TT += p.P4();
    }
    
    double const partonMtt = p4TT.M();
    
    
    // Apply smearing
    double const smearedMtt = rGen.Gaus(partonMtt, partonMtt * resolution);
    
    
    // Fill nominal histogram and histograms with mtt scale variations
    double const nominalWeight = reader->GetWeight();
    histNominal->Fill(smearedMtt, nominalWeight);
    histScaleUp->Fill(smearedMtt * (1 + scaleVariation), nominalWeight);
    histScaleDown->Fill(smearedMtt * (1 - scaleVariation), nominalWeight);
    
    
    auto const &lheWeights = reader->GetLHEWeights();
    
    // Book histograms with alternative weights if not done yet
    if (not histAltWeightsBooked)
    {
        for (auto const &weight: lheWeights)
            histAltWeights.emplace_back(
                processor->Book<TH1D>("", ("AltWeight_ID" + std::to_string(weight.ID)).c_str(), "",
                    binning.size() - 1, binning.data())
            );
        
        histAltWeightsBooked = true;
    }
    
    
    // Fill histograms with alternative weights
    for (unsigned i = 0; i < lheWeights.size(); ++i)
        histAltWeights[i]->Fill(smearedMtt, lheWeights[i].Weight);
    
    
    return true;
}
