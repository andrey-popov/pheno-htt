#include <DelphesReaderBase.hpp>


DelphesReaderBase::DelphesReaderBase(unsigned readOptions):
    jetPtThreshold(20.), jetEtaThreshold(2.4)
{
    readLHEParticles = ((readOptions & LHE_PARTICLES) == LHE_PARTICLES);
}
