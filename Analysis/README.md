# Analysis

In this section phenomenological analysis is performed.
SM tt background is described using MC simulation, as documented in section &ldquo;[Generation](../Generation)&rdquo;, while signal distributions are constructed on the fly.

Dependencies:

 * [ROOT](root.cern.ch) 6.14/04.
 * [LHAPDF](https://lhapdf.hepforge.org/) 6.1.6 with Python bindings. Download PDF set `PDF4LHC15_nlo_30_pdfas` (LHAPDF ID 90400).
 * Python 3.5.
 * [NumPy](http://numpy.org) 1.13.1, [SciPy](https://scipy.org/scipylib/index.html) 0.18.1, [Matplotlib](https://matplotlib.org) 2.0.2.


## SM tt

Produce m<sub>tt</sub> histograms for SM tt, both nominal and systematic variations, by running

```sh
./smHists.sh
./buildTemplates.py
rm -r hists
```

Internally, this executes program [mtt-hists](../EventProcessing/prog/mtt-hists.cpp) for SM tt ROOT files produced in section &ldquo;[Generation/Showering](../Generation/Showering)&rdquo;.
A default smearing of 20% is applied to mimic reconstruction effects; its value can be configured if needed.
Script [`buildTemplates.py`](buildTemplates.py) combines all histograms in a single file `ttbar.root`.
It also applies a k-factor of 2 and efficiency of additional event selection.
All templates are normalized to an integrated luminosity of 1&thinsp;pb<sup>-1</sup>.
Constructed systematic variations can be plotted using script [`plotVariations.py`](plotVariations.py).

Systematic variations that are not described by per-event weights suffer from statistical fluctuations, which can pose a problem for the subsequent statistical analysis.
To protect against this, they are smoothed with a version of the [LOWESS algorithm](https://en.wikipedia.org/wiki/Local_regression).
Run with

```sh
./smoothTemplates.py
mv ttbar_smooth.root ttbar_res20.root
rm ttbar.root
```

Script `smoothTemplates.py` also plots the resulting variations.
Produced file `ttbar_res20.root` with nominal m<sub>tt</sub> in SM tt and its systematic variations is [stored](ttbar_res20.root) in the repository for convenience, as well as the same file produced with 10% resolution, [`ttbar_res10.root`](ttbar_res10.root).


## Signal

Signal templates are constructed from analytical parton-level cross sections, which are provided (both for the resonant part and the interference) for various signal models. For example, module [`twohdm.py`](twohdm.py) implements a generic 2HDM, which is then specialized in [`hmssm.py`](hmssm.py) to describe hMSSM. Typically, k-factors are applied to account for higher-order corrections in the gg&Phi; form-factor.

In order to accelerate the convolution of the parton-level cross sections with the PDF, one of the involved integrals over Bjorken&nbsp;x of the incoming gluons is precomputed with

```sh
./sHatPDF.py
```

storing the result in file `sHatPDF.npy`, which is also [included](sHatPDF.npy) in the repository for convenience.

Module [`spectrum.py`](spectrum.py) provides a class to compute differential cross section in reconstructed m<sub>tt</sub> starting from the parton-level cross section.
It applies the branching ratio for the targeted decays (set to 8/27, which corresponds to &ell;&thinsp;+&thinsp;jets), performs the convolution with PDF, applies efficiency of the event selection (a hard-coded function of parton-level m<sub>tt</sub>), and performs the smearing to account for the resolution effects.

The full setup can be tested with script [`plotSmearedSignal.py`](plotSmearedSignal.py), which constructs and plots m<sub>tt</sub> distributions in hMSSM for different values of the resolution.


## Parameter scans

Parameter scans for various signal models are performed with scripts `scan_*.py`. Each one can be run for different experimental benchmarks, which are defined by the integrated luminosity and m<sub>tt</sub> resolution (which also determines what file with background templates should be used). Here is an example for hMSSM:

```sh
./scan_hMSSM.py --bkg ttbar_res20.root --resolution 0.2 --lumi 150 -o fig/hMSSM_res20_150ifb.pdf
```

Due to a large number of points, this command takes about an hour.
For each probed point on the (m<sub>A</sub>, tan&thinsp;&beta;) plane, the script computes the expected significance and the CL<sub>s</sub> value for an upper limit.
