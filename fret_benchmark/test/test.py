#!/usr/bin/env python

import unittest
import math
import collections
import sys
import os
import re
import subprocess

scripts = 'DATA/SCRIPTS'
ref = 'DATA/REFERENCE_RESULTS'

Data = collections.namedtuple('Data', 'mean stddev')

def get_data(fname):
    r = re.compile('(\w+)\s+([\d\.]+) \( ([\d\.]+) \)')
    line = open(fname).readline()
    data = {}
    for m in r.finditer(line):
        # Note: convert standard error to stddev by multiplying by sqrt(ntests)
        data[m.group(1)] = Data(mean=float(m.group(2)),
                                stddev=float(m.group(3)) * math.sqrt(20))
    return data

class Tests(unittest.TestCase):

    def test_linker(self):
        """Test with linker flexibility"""
        self.run_test("linker")

    def test_nolinker(self):
        """Test without linker flexibility"""
        self.run_test("nolinker")

    def run_test(self, testtype):
        # Just run one test
        p = subprocess.check_call(['./do_test.sh', '1'],
                                  cwd='%s/BENCHMARK_%s' % \
                                  (scripts, testtype.upper()))
        # Analysis
        p = subprocess.check_call(['%s/ANALYSIS/do_analysis.py' % scripts,
                                   'RESULTS_%s' % testtype.upper(),
                                   '%s.dat' % testtype])
        p = subprocess.check_call(['%s/ANALYSIS/do_analysis_per_complex.py' \
                                   % scripts,
                                   'RESULTS_%s' % testtype.upper(),
                                   '%s_per_complex.dat' % testtype])
        # The predicted mean (= the value from the single test) should lie
        # within a few standard deviations of that in the reference output
        data = get_data('%s_per_complex.dat' % testtype)
        ref_data = get_data('%s/%s_per_complex.dat' % (ref, testtype))
        for key in ref_data.keys():
            self.assertNearMean(key, data[key].mean, ref_data[key].mean,
                                ref_data[key].stddev, 5)

    def assertNearMean(self, key, mean, ref_mean, ref_stddev, mult):
        self.assert_(abs(mean - ref_mean) < ref_stddev * mult,
                     "Predicted mean for %s (%.2f) is not within %d "
                     "standard deviations (%d * %.2f) of the expected "
                     "mean (%.2f)" % (key, mean, mult, mult,
                                      ref_stddev, ref_mean))

if __name__ == '__main__':
    # Always run from top-level directory
    os.chdir(os.path.join(os.path.dirname(sys.argv[0]), '..'))
    unittest.main()
