# Phenomenological study of &Phi;&rarr;tt

[![DOI](https://zenodo.org/badge/117536495.svg)](https://zenodo.org/badge/latestdoi/117536495)

This repository provides supporting material for a phenomenological study of the production of a heavy Higgs boson decaying into a pair of top quarks: A.&nbsp;Djouadi, J.&nbsp;Ellis, A.&nbsp;Popov, J.&nbsp;Quevillon, *Interference effects in tt production at the LHC as a window on new physics* [JHEP 03 (2019) 119](https://doi.org/10.1007/JHEP03(2019)119) [[arXiv:1901.03417](https://arxiv.org/abs/1901.03417)]. The code used to perform the scans over parameters of various signal models is included here. It is organized into three blocks:

 * [Generation](Generation): Production of LHE files for several signal benchmarks and SM tt background, showering and hadronization.
 * [Event processing](EventProcessing): Programs to work with produced files with simulated events.
 * [Analysis](Analysis): Scripts to perform parameter scans for various signal models and to produce inputs required for them.

Details for each block are provided in dedicated README files. Readers interested in reproducing the parameter scans can proceed directly to [this section](Analysis#parameter-scans). All input files needed for the scans are included in the repository.
