# FRET Benchmark

This directory contains all the scripts to run the benchmark of the FRETR
Bayesian restraint, and the reference results.  

## Folder structure 

- `DATA`:
  - `GMM`:  Gaussian Mixture Models needed by model for linker flexibility
  - `PDBs`: Pdb with complex structure, one per chain
  - `REFERENCE_RESULTS`: reference results of benchmark
  - `SCRIPTS`: scripts to run and analyze the benchmark
     - `ANALYSIS`: scripts to analyze the results
     - `BENCHMARK_LINKER`: scripts to run the benchmark with model for linker flexibility
     - `BENCHMARK_NOLINKER`: scripts to run the benchmark without model for linker flexibility
      
- `RESULTS_LINKER`   : directory where the results of benchmark with model for linker flexibility will be stored
- `RESULTS_NOLINKER` : directory where the results of benchmark without model for linker flexibility will be stored

## How to run the benchmark 

0) This is the `ROOT` directory.

1) PRELIMINARY. Compile IMP (any recent version should work - see below for the last known good IMP version). IMP git version develop-c47408c was used for the publication.

2) EXECUTION. Go to `DATA/SCRIPTS/BENCHMARK_LINKER` or `DATA/SCRIPTS/BENCHMARK_NOLINKER` where you will find
   the scripts to run the benchmark with or without using a model for linker flexibility, respectively. Then: 

   - To run a specific test, run:

      `./do_test.sh ID`

     where `ID` is an integer from 1 to 1280.

   - To run all tests, just use a bash loop:

      `for((ID=1;ID<=1280;ID++)); do ./do_test.sh $ID; done`

   - A sample SGE script to run all the tests on a cluster is also provided, as `job.sh`.

3) ANALYSIS. From the `ROOT` directory, do:

   - `DATA/SCRIPTS/ANALYSIS/do_analysis.py             RESULTS_LINKER    linker.dat`    
   - `DATA/SCRIPTS/ANALYSIS/do_analysis.py             RESULTS_NOLINKER  nolinker.dat`
   - `DATA/SCRIPTS/ANALYSIS/do_analysis_per_complex.py RESULTS_LINKER    linker_per_complex.dat`
   - `DATA/SCRIPTS/ANALYSIS/do_analysis_per_complex.py RESULTS_NOLINKER  nolinker_per_complex.dat`

   Statistics will be written to the `.dat` files, which should be compared to the results in `DATA/REFERENCE_RESULTS` and to those
   published in Bonomi et al. "Protein complex structures from Bayesian modeling of in vivo FRET data", as follows.
   - `linker.dat`             -> Table S2
   - `nolinker.dat`           -> Table S3


_Author(s)_: Max Bonomi

_Version_: 1.0


_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Last known good IMP version_: [![build info](https://salilab.org/imp/systems/?sysstat=3)](http://salilab.org/imp/systems/)

_Testable_: Yes.

_Parallelizeable_: Yes

_Publications_:
 - Bonomi et al. "Protein complex structures from Bayesian modeling of in vivo FRET data", Submitted 2013. 
