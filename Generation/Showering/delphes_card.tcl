#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  NeutrinoFilter
  GenJetFinder
  GenMissingET

  TreeWriter
}



#####################
# Neutrino Filter
#####################

module PdgCodeFilter NeutrinoFilter {

  set InputArray Delphes/stableParticles
  set OutputArray filteredParticles

  set PTMin 0.0

  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}

}

#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 20.0
}

#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}


##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
#  add Branch Delphes/allParticles Particle GenParticle
  add Branch Delphes/allParticlesLHEF ParticleLHEF GenParticle

  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenMissingET/momentum GenMissingET MissingET
}
