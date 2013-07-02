# FRET Benchmark

This directory contains all the scripts to run the benchmark of the FRETR
Bayesian restraint, and the reference results.  
The benchmark should be run with IMP git version develop-c47408c.

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

1) PRELIMINARY. Compile IMP git version develop-c47408c

2) EXECUTION. Go to `DATA/SCRIPTS/BENCHMARK_LINKER` or `DATA/SCRIPTS/BENCHMARK_NOLINKER` where you will find
   the scripts to run the benchmark with or without using a model for linker flexibility, respectively. Then: 

   - To run a specific test, first open `do_test.sh` and set the `IMP` and `ROOT` directories.
     Please, do not change the name of the other directories. Then:

      `bash do_test.sh ID`

     where `ID` is an integer from 1 to 1280.

   - To run all tests, just use a bash loop:

      `for((ID=1;ID<=1280;ID++)); do bash do_test.sh $ID; done`

   - To run all tests on the QB3 cluster, first open `job.sh` and set the `IMP` and `ROOT` directories, as well as
     the name of the directories where standard output and error will be written.
     Please, do not change the name of the other directories. Then:

      `qsub job.sh`

3) ANALYSIS. From the `ROOT` directory, do:

   - `python DATA/SCRIPTS/ANALYSIS/do_analysis.py             RESULTS_LINKER    linker.dat`    
   - `python DATA/SCRIPTS/ANALYSIS/do_analysis.py             RESULTS_NOLINKER  nolinker.dat`
   - `python DATA/SCRIPTS/ANALYSIS/do_analysis_per_complex.py RESULTS_LINKER    linker_per_complex.dat`
   - `python DATA/SCRIPTS/ANALYSIS/do_analysis_per_complex.py RESULTS_NOLINKER  nolinker_per_complex.dat`

   Statistics will be written to the `.dat` files, which should be compared to the results in `DATA/REFERENCE_RESULTS` and to those
   published in Bonomi et al. "Protein complex structures from Bayesian modeling of in vivo FRET data", as follows.
   - `linker.dat`             -> Table S1
   - `nolinker.dat`           -> Table S2
   - `linker_per_complex.dat` -> Table S3
   - `linker_per_complex.dat` -> Table S4


_Author(s)_: Max Bonomi

_Version_: 1.0


_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Testable_: Yes.

_Parallelizeable_: Yes

_Publications_:
 - Bonomi et al. "Protein complex structures from Bayesian modeling of in vivo FRET data", Submitted 2013. 
