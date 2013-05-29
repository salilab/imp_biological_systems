#!/usr/bin/env python

from modeller import *

env = environ()
aln = alignment(env)
#directories for input atom files
env.io.atom_files_directory = ['.', 'data/templates']
mdl = model(env, file='1iok', model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='1iok', atom_files='1iok.pdb')
aln.append(file='data/sequences/groel_ecoli.ali', align_codes='P0A6F5')
aln.align2d(max_gap_length=50)
aln.write(file='groel-1iokA.ali', alignment_format='PIR')
