# PDE6

These scripts demonstrate the use of IMP,
[MODELLER](http://salilab.org/modeller/) and
[Chimera](http://www.cgl.ucsf.edu/chimera/) in the modeling of the
phosphodiesterase (PDE6).

## Homology modeling

First, MODELLER is used to generate initial model for the structure using multople templates (PDB codes 1tbf, 3bjc, 3dba, 3ibj). Cross links, symmetry and secondary structure restraints are used in modelling:
 `cd comparative_modelling; ./model_mult.py`

The resulting models do not fit well into 3D EM density map (20A resolution).

Therefore, IMP is used to fit the complex into the electron microscopy density map.



## Integrative modeling with IMP

It uses IMP svn Revision: 16273

with the following patch:

```
Index: modules/em/src/EnvelopeFitRestraint.cpp
===================================================================
--- modules/em/src/EnvelopeFitRestraint.cpp	(revision 16273)
+++ modules/em/src/EnvelopeFitRestraint.cpp	(working copy)
@@ -45,9 +45,9 @@
   double best_score = -std::numeric_limits<double>::max();
   for(unsigned int j=0; j<map_transforms.size(); j++) {
     //std::cerr << "Scoring " << map_transforms[j] << std::endl;
-    if(!envelope_score_.is_penetrating(coordinates,
-                                       map_transforms[j],
-                                       penetration_thr)) {
+    //if(!envelope_score_.is_penetrating(coordinates,
+    //                                   map_transforms[j],
+    //                                   penetration_thr)) {
       //std::cerr << "  not penetrating " << map_transforms[j] << std::endl;
       double score = envelope_score_.score(coordinates, map_transforms[j]);
       //std::cerr << "  score = " << score << std::endl;
@@ -56,7 +56,7 @@
         best_trans = map_transforms[j];
         best_found = true;
       }
-    }
+    //}
   }

   if(best_found)
```
Usage and content of the directory `integrative_modeling`

1)  test the python script:
`./run_modeling.py test`
no error = all tests passed

2) run sampling:
`./run_modeling.py`
all output data will be stored in output/

3) analyse the results:

3.1) compile the clustering algorithm in `bin/`
`gfortran cluster.f u3best.f -o cluster.x`

3.2) run the analysis script:
`bin/get_frames.sh`

100 best scoring frames will be stored in `best_pdb/`
clustering data will be stored in `clustering/`

## Refinement

Finally the final models are rebuilt and refined with modeller:
`model_refinement/cluster1/model-single.py`;


_Author(s)_: Riccardo Pellarin, Dina Schneidman

_Version_: 1.0


_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Testable_: Yes.

_Parallelizeable_: No

_Publications_:
 - in preparation
