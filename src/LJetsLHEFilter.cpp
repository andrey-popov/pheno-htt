#include <LJetsLHEFilter.hpp>

#include <DelphesReaderBase.hpp>

#include <cmath>


LJetsLHEFilter::LJetsLHEFilter(DelphesReaderBase const *reader_):
    reader(reader_)
{}


bool LJetsLHEFilter::ProcessEvent()
{
    unsigned nE = 0, nMu = 0, nTau = 0;
    
    for (auto const &p: reader->GetLHEParticles())
    {
        int const absPdgID = std::abs(p.PID);
        
        if (absPdgID == 11)
            ++nE;
        else if (absPdgID == 13)
            ++nMu;
        else if (absPdgID == 15)
            ++nTau;
    }
    
    
    return ((nE + nMu) == 1 and nTau == 0);
}
