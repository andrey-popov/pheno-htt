# Generation

In this section SM&nbsp;tt and &Phi;&rarr;tt events are generated.


## LHE tt

Produce LHE files for SM&nbsp;tt events using [MadGraph5_aMC@NLO](https://launchpad.net/mg5amcnlo) program.

Events are generated at leading order, with no additional jets in the final state.
Renomalization and factorization scales are set dynamically to m<sub>tt</sub>&nbsp;/&nbsp;2.
Use NLO [PDF4LHC15](https://arxiv.org/abs/1510.03865) PDF set.
Top quarks are fully decayed within MadGraph, and only semileptonic final states containing a muon or an electron are kept.
Events in which charged leptons have p<sub>T</sub>&nbsp;&gt;&nbsp;30&nbsp;GeV are selected.

Dependencies:

 * [LHAPDF](https://lhapdf.hepforge.org/) 6.1.6. Download PDF set `PDF4LHC15_nlo_30_pdfas` (LHAPDF ID 90400).
 * [MadGraph5_aMC@NLO](https://launchpad.net/mg5amcnlo) 2.6.0. Madgraph executables must be available in the `$PATH`. Copy provided file [`setscales.f`](LHE-tt/setscales.f) into `MG5_aMC_v2_6_0/Template/LO/SubProcesses`. This file provides a user-defined scale set as described above. MadGraph must be configured to use LHAPDF.
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

For each mass 5M events will be produced, with 100k events per script.
These are then launched with

```sh
./runMadGraph.sh configs/*
```

Produced LHE events contain alternative weights for systematic variations in the scales and PDF.

For convenience an example MadGraph banner is [provided](LHE-tt/ttbar_banner.txt) in the repository.


## LHE signal

Produce LHE files for signal events (both the resonant part and the interference).
The setup and the event selection are the same as for the SM&nbsp;tt above.

Subdirectory [`HeavyHiggs`](LHE-signal/HeavyHiggs) contains the MadGraph model for the signal process.
It is based on the SM and adds CP-odd and even Higgs bosons.
Their masses, widths, and couplings to top quarks are adjustable parameters of the model.
An effective coupling to gluons is implemented with a top quark loop following [Spira et al.](https://arxiv.org/abs/hep-ph/9504378)
The model was originally developed by Sébastien Brochet and Stéphane Perriès, and it is copied from [here](https://github.com/IPNL-CMS/HTTMadgraphDocumentation/tree/19854e06bae1e861637b20904f6f6db21f8a42fc/model/Massive_Higgs_UFO).
To install it, the directory [`HeavyHiggs`](LHE-signal/HeavyHiggs) needs to be copied into `MG5_aMC_v2_6_0/models`.

For the resonant part, events are generated in a similar way to SM&nbsp;tt:

```sh
mg5_aMC A-res-setup.script
./writeConfigs.py A A-res -o configs/A-res
./runMadGraph.sh configs/A-res/*
```

and similarly for the CP-even state H. It will produce 1M events (with 100k events per LHE file) for each of five mass points specified in script [`writeConfigs.py`](LHE-signal/writeConfigs.py), with the width set to 10%.

On the other hand, generation of events for the interference part with decayed top quarks requires working around a [limitation](https://bugs.launchpad.net/mg5amcnlo/+bug/1511378) of the current version of MadGraph5_aMC@NLO.
As suggested by Sébastien Wertz, the matrix element for the full SM+BSM process is constructed and then the corresponding Fortran code is modifed to keep only terms relevant to the interferences.
This is done with script [`patchMEInt.py`](LHE-signal/patchMEInt.py).
Here is a summary of commands for the CP-odd state:

```sh
mg5_aMC A-int-setup.script
./patchMEInt.py A-int 85
./writeConfigs.py A A-int -o configs/A-int
./runMadGraph.sh configs/A-int/*
```

In the above, 85 is the index of the [Att coupling](LHE-signal/HeavyHiggs/couplings.py#L348).
For the CP-even state the index is [86](LHE-signal/HeavyHiggs/couplings.py#L352).


## Showering

Perform showering and hadronization of produced LHE files.
Here the SM&nbsp;tt is used as the example, while signal files can be processed in the same way.

Showering and hadronization are done with [Pythia 8](http://home.thep.lu.se/~torbjorn/Pythia.html) program, using the default tune Monash 2013.
For convenience, it is run inside [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes) framework, which also clusters stable particles into jets (using the [FastJet](http://fastjet.fr) package) and saves events in a binary format.
A [minimalistic configuration](Showering/delphes_card.tcl) for Delphes is used, which only saves in the output file LHE-level particles and jets.

Dependencies:

 * [ROOT](root.cern.ch) 6.10/04.
 * [Pythia](http://home.thep.lu.se/~torbjorn/Pythia.html) 8.230.
 * [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes) 3.4.1 built with support of Pythia.
 * Python 3.5.
 
Process LHE files with the following commands:

 ```sh
 ./runDelphes.py ../LHE/SM-tt/Events/ttbar_* -o SM-tt
 ./runDelphes.py ../LHE/SM-tt/Events/ttbar_[0-9]* --pythia-config pythiaConfig_FSR-up.cmnd -o SM-tt_FSR-up
 ./runDelphes.py ../LHE/SM-tt/Events/ttbar_[0-9]* --pythia-config pythiaConfig_FSR-down.cmnd -o SM-tt_FSR-down
 ```
 
Here [runDelphes.py](Showering/runDelphes.py) is a convenience script that runs multiple copies of the `DelphesPythia8` executable in parallel and takes care of unpacking compressed LHE files and dealing with temporary files.
LHE files with nominal m<sub>t</sub> are processed three times with the renormalization scale used in showering set to its nominal value or scaled by factors 0.5 and 2.
After the generation is done, log files can be removed and Delphes files from directories `SM-tt_FSR-*` can be renamed and moved into `SM-tt` by running script `cleanup.sh`.
