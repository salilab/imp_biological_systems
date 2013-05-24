#!/usr/bin/env python

import subprocess

cmd = ['multifit.py', 'segment', "data/em_maps/groel-11.5A.mrc",
       "14", "0.7", "output/groel_segments_center.pdb",
       "--seg", "output/groel_subunit", "--apix", "2.7",
       "--x", "-135", "--y", "-135", "--z", "-135"]
print "======running command:"
print " ".join(cmd)
subprocess.call(cmd)
