# Generation

In this section SM tt events are generated.


## LHE

Produce LHE files for SM tt events using [MadGraph5_aMC@NLO](https://launchpad.net/mg5amcnlo) program.

Events are generated at leading order, with no additional jets in the final state.
Renomalization and factorization scales are set dynamically to m<sub>tt</sub>&nbsp;/&nbsp;2.
Use NLO [PDF4LHC15](https://arxiv.org/abs/1510.03865) PDF set.
Top quarks are fully decayed within MadGraph, and only semileptonic final states containing a muon or an electron are kept.
Events in which charged leptons have p<sub>T</sub>&nbsp;&gt;&nbsp;30&nbsp;GeV are selected.

Dependencies:
 * [LHAPDF](https://lhapdf.hepforge.org/) 6.1.6. Download PDF set `PDF4LHC15_nlo_30_pdfas` (LHAPDF ID 90400).
 * [MadGraph5_aMC@NLO](https://launchpad.net/mg5amcnlo) 2.6.0. Madgraph executables must be available in the `$PATH`. Copy provided file [`setscales.f`](LHE/setscales.f) into `MG5_aMC_v2_6_0/Template/LO/SubProcesses`. This file provides a user-defined scale set as described above. MadGraph must be configured to use LHAPDF.
 * Python 2.7 (required by MadGraph).

To generate events, first create a MadEvent directory for the process:
```sh
mg5_aMC SM-tt-setup.script
```
Create scripts that will, starting from the MadEvent directory, generate SM tt events with three different values of m<sub>t</sub> (nominal and a systematic variation):
```sh
./writeConfigs.py SM-tt --mt 173 -o configs/ttbar
./writeConfigs.py SM-tt --mt 173.5 -o configs/ttbar_mt-up
./writeConfigs.py SM-tt --mt 172.5 -o configs/ttbar_mt-down
```
For each mass 5M events will be produced, with 100k events per script. The generation is then launched with
```sh
./runMadGraph.sh configs/*
```

Produced LHE events contain alternative weights for systematic variations in the scales and PDF.

For convenience an example MadGraph banner is [provided](LHE/ttbar_banner.txt) in the repository.
