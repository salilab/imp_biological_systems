# Nup84

The nuclear pore complex (NPC) is a multiprotein assembly that serves as the sole mediator
of nucleocytoplasmic exchange in eukaryotic cells. Here, we use an integrative approach to determine
the structure of an essential component of the yeast NPC, the ~600 kDa heptameric Nup84 complex. The configuration of the subunit structures was determined by satisfaction of
spatial restraints derived from a diverse set of negative stain EM and protein domain mapping data.
Phenotypic data was mapped onto the complex, allowing us to find regions that stabilize the NPC's
interaction with the NE membrane and connect the complex to the rest of the NPC. We suggest a
scenario for the evolution of the Nup84 complex through a series of gene duplication and loss events.
This work demonstrates that integrative approaches based on low resolution data can generate
functionally informative structures at intermediate resolution.

## How to run

The commands given here are for running in the Nup84 directory.

## Input files
- `input_models.txt`: a list of input models and their pdb
file names, plus a which subsequence of which protein each model
represents
- `input/structures`: the 14 input model pdb files
- `nup84_class_average.jpg`: contains a 2d class average of
particles from negative stain em located in input/em
- `coarse_restraints_set1`,
`coarse_restraints_set`, `coarse_restraints_set3`, `refined_restraints`
located in `input/restraints` are the restraint data input files.

## Generating models

The following steps generate one model. They should be repeated many times to
generate an ensemble of models. When %1% occurs in an argument, it should be
replaced by a number or something similar to differentiate the various runs.

First, generate a coarse structural solution
`../../tools/imppy.sh python coarse_stage1.py input_models.txt input/restraints coarse_structures/%1%.pdb coarse_structures/%1%_info.txt`

Second, refine the model
`../../tools/imppy.sh python refine_stage2.py input_models.txt
input/restraints coarse_structures/%1%.pdb coarse_structures/%1%_info.pdb refined_structures/%1%.pdb refined_structures/%1%_info.txt`

Once you have set of refined models, generate em2d scores for all of them
`../../tools/imppy.sh python score_em2d.py nup84_class_average.jpg refined_structures refined_structures_em2d_scores.txt`

Next, filter out all except the top 10% of the structures based  on the em2d scores
`../../tools/imppy.sh python filter_stage3.py refined_structures/ filtered_structures/  refined_structures_em2d_scores.txt  10`


Next, proceed with the final refinement of each of the filtered structures
`../../tools/imppy.sh  python refine_stage4.py input_models.txt input/restraints nup84_class_average.jpg filtered_structures/%1%.pdb final_structures/%1%.pdb final_structures/%1%_info.txt`

## Analysis

Finally, generate 7 localization density maps from the final set of
structures (one per protein). Any of the `info` files from the last
step will do.
`../../tools/imppy.sh python generate_dmap.py final_structures/ final_structures/%1%_info.txt localization_maps/`


## Info

_Author(s)_: Jeremy Phillips

_Version_: 1.0


_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Testable_: No.

_Parallelizeable_: Yes.

_Last known good \imp version_: None

_Publications_:
 - J. Fernandez-Martinez, J. Phillips, M. Sekedat, R. Diaz-Avalos, J. Velazquez-Muriel, J. Franke, R. Williams, D. Stokes, B. Chait, A. Sali, M. Rout, Structure-function Map for a Heptameric Component of the Nuclear Pore Complex, 2011.
