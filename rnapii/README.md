# Rnapii

\note Currently, generating the structure requires \c IMP.multifit2. It can be obtained from the author (Keren). We are working on moving the needed functionallity into \imp propper. Your patience is appreciated.

## How to run

All command here are given as if run from the \c biological_systems/rnapii directory in the build directory

## Input

 - `rnapii.subunit.list.txt`: The subunits list. Each line describes a subunit and contains: `name pdb_file [global/local]_fit`

 - `data/emd_1283.mrc`: a mrc file of the assembly

 - `data/*.pdb`: pdb files of the subunits, derived from the yeast homologues.

 - `rnapii.subunit.list.txt`: the list of subunits and corresponding pdb files

## Generating the parameters file

`../../tools/imppy.sh python ../../modules/multifit2/bin/generate_assembly_input.py -i rnapii.asmb.input -- rnapii data/rnapii.subunit.list.txt 30 data/emd_1283.mrc 20 2.5 0.2 60 60 60`

This creates the file `rnapii.asmb.input` which can be used to run multifit.

## Running multifit
First, create all the molecular surfaces to be used when computing excluded volumes.
`../../tools/imppy.sh python ../../modules/multifit2/bin/create_all_surfaces.py rnapii.asmb.input`

Creates `X.ms` for each input pdb.

Second, generate the assembly anchor graph:
`../../tools/imppy.sh python ../../modules/multifit2/bin/generate_fine_assembly_anchor_graph.py rnapii.asmb.input`

This creates:
- `rnapii.asmb.anchors.pdb`  : The graph in pdb format
- `rnapii.asmb.anchors.txt`  : The graph in txt format
- `rnapii.asmb.anchors.cmm`  : The graph in cmm format

Then generate the fits
`../../tools/imppy.sh python $IMP/modules/multifit2/bin/run_fitting_fft.py -p rnapii.multifit.param  rnapii.asmb.data -c 6 -z -n 10000`
The fits are put in files named \c X_fitting.txt.

Then generate the indexes
`../../tools/imppy.sh python $IMP/modules/multifit2/bin/generate_indexes_from_fitting_solutions.py rnapii rnapii.asmb.input 1000 -f data/subunit.indexes`
The indexes are put into files name `X.fit.indexes.txt`.

Next we create the automated proteomics file
`../../tools/imppy.sh python ../../modules/multifit2/bin/create_auto_proteomics_file.py rnapii.asmb.input rnapii_em_coarse_anchors_FINE.txt rnapii.ev_proteomics.input`

This produces the file `rnapii.ev_proteomics.input` for assembling the components using only excluded volume.

Finally we create the assembly fitting structures.
`../../tools/imppy.sh python $IMP/modules/multifit2/bin/align_proteomics_em_atomic_plan.py -m 30 rnapii.asmb.input rnapii.ev_proteomics.input rnapii.indexes.mapping.input rnapii.alignment.param rnapii.docking.param rnapii.combinations.output rnapii.scores ignored`

If you want to use proteomics restraints taken from [BioGrid](http://thebiogrid.org), you can modify the \c rnapii.ev_proteomics.input file, or use the provided `rnapii.biogrid_proteomics.input` instead.

The output files
- `rnapii.combinations.output` gives the indexes into `X_fitting.txt` of each solution.
- `rnapii.scores` gives the fitting score of each solution.

## Analyzing and displaying the results

To convert the output files to pdbs do
`../../tools/imppy.sh python ../../modules/multifit2/bin/write_ensemble_models.py 1z5s.asmb.input combinations.output asmb.mdl`

And to print more detailed information about the various excluded volume,
docking and proteomics scores run:
`../../tools/imppy.sh ../../modules/multifit2/bin/score_ensemble_modles.py rnapii.asmb.input rnapii.proteomics.input rnapii.indexes.mapping.input rnapii.alignment.param rnapii.docking.param rnapii.combinations.output rnapii_detailed.scores`


## Info

_Author(s)_: Keren Lasker

_Version_: 1.0


_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Testable_: No.

_Parallelizeable_: No.

_Last known good IMP version_: None

_Publications_:
 - Keren Lasker, Daniel Russel, Jeremy Phillips, Haim Wolfson, Andrej Sali, Determining architectures of macromolecular assemblies by aligning interaction networks to electon microscopy density maps, 2011.
