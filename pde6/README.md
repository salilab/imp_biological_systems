# PDE6

These scripts demonstrate the use of IMP,
[MODELLER](http://salilab.org/modeller/) and
[Chimera](http://www.cgl.ucsf.edu/chimera/) in the modeling of the
phosphodiesterase (PDE6).

First, MODELLER is used to generate initial model for the structure using multople templates (PDB codes 1tbf, 3bjc, 3dba, 3ibj). Cross links, symmetry and secondary structure restraints are used in modelling:
 `comparative_modelling/model_mult.py`

The resulting models do not fit well into 3D EM density map (20A resolution).

Therefore, IMP is used to fit the complex into the electron microscopy density map.

...

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
