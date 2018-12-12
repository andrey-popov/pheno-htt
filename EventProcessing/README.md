# Event processing

This directory contains source code for several programs to process ROOT files produced in section &ldquo;[Generation/Showering](../Generation/Showering)&rdquo;. Only one of the programs, [mtt-hists](prog/mtt-hists.cpp), was used for results included in the paper. It is run by scripts in section &ldquo;[Generation/Showering](../Analysis)&rdquo;.

Dependencies:

 * CMake 2.8.12.
 * C++ compiler supporting C++14 standard.
 * [ROOT](root.cern.ch) 6.14/04.
 * [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes) 3.4.1.
