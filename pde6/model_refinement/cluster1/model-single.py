#!/usr/bin/env python

import sys
from modeller import *
from modeller.automodel import *

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms

        # secondary structure based on Quick2D predictions

        # secondary structure - chain A
        rsr.add(secondary_structure.alpha(self.residue_range('6:A', '13:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('18:A', '24:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('31:A', '39:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('42:A', '69:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('288:A', '296:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('372:A', '377:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('485:A', '491:A')))

        # secondary structure - chain B
        rsr.add(secondary_structure.alpha(self.residue_range('6:B', '12:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('19:B', '23:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('29:B', '37:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('42:B', '68:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('285:B', '293:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('370:B', '375:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('483:B', '490:B')))

        # cross links
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:21:A'], at['CA:31:A']), mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:76:A'], at['CA:83:A']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:76:A'], at['CA:84:A']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:534:A'], at['CA:581:A']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:765:A'], at['CA:827:A']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:819:A'], at['CA:827:A']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:677:A'], at['CA:683:A']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:827:A'], at['CA:845:A']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:786:A'], at['CA:807:A']),  mean=12.0, stdev=0.5))

        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:784:B'], at['CA:805:B']), mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:440:B'], at['CA:392:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:490:B'], at['CA:532:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:532:B'], at['CA:579:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:532:B'], at['CA:577:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:675:B'], at['CA:681:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:687:B'], at['CA:784:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:763:B'], at['CA:817:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:817:B'], at['CA:832:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:817:B'], at['CA:826:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:817:B'], at['CA:839:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:819:B'], at['CA:832:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:826:B'], at['CA:833:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:832:B'], at['CA:839:B']),  mean=12.0, stdev=0.5))

        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:21:A'], at['CA:26:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:21:A'], at['CA:44:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:31:A'], at['CA:78:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:31:A'], at['CA:26:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:247:A'], at['CA:391:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:447:A'], at['CA:440:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:447:A'], at['CA:618:B']),  mean=12.0, stdev=0.5))

        # Constrain the A and B chains to be identical (but only restrain
        # the C-alpha atoms, to reduce the number of interatomic distances
        # that need to be calculated):
        s1 = selection(a.residue_range('12:A', '183:A')).only_atom_types('CA') | selection(a.residue_range('225:A', '481:A')).only_atom_types('CA') | selection(a.residue_range('508:A', '624:A')).only_atom_types('CA') | selection(a.residue_range('629:A', '802:A')).only_atom_types('CA')
        s2 = selection(a.residue_range('10:B', '181:B')).only_atom_types('CA') | selection(a.residue_range('223:B', '479:B')).only_atom_types('CA') | selection(a.residue_range('506:B', '622:B')).only_atom_types('CA') | selection(a.residue_range('627:B', '800:B')).only_atom_types('CA')

        self.restraints.symmetry.append(symmetry(s1, s2, 1.0))
#     def user_after_single_model(self):
#         # Report on symmetry violations greater than 1A after building
#         # each model:
#         self.restraints.symmetry.report(1.0)
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['A', 'B'],
                             renumber_residues=[1, 1])


env = environ()
a = MyModel(env, alnfile='pde6-pde6.ali',
            knowns='model.pdb', sequence='pde6',
            assess_methods=(assess.DOPE, assess.GA341))
#a.very_fast()                       # prepare for extremely fast optimization
a.starting_model = 1
a.ending_model = 5

# If testing is requested, build just one model:
if sys.argv[-1] == '--run_quick_test':
    a.ending_model = 1

a.make()
