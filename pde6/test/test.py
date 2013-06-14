#!/usr/bin/env python

import unittest
import sys
import os
import subprocess

class Tests(unittest.TestCase):

    def test_comparative_modeling(self):
        """Test the comparative modeling script"""
        # Make sure the script runs without errors
        p = subprocess.check_call(['./model_mult.py', '--run_quick_test'],
                                  cwd='comparative_modelling')

    def test_integrative_modeling(self):
        """Test the integrative modeling script"""
        # Make sure the script runs without errors
        p = subprocess.check_call(['./run_modeling.py', 'test'],
                                  cwd='integrative_modeling')

if __name__ == '__main__':
    # Always run from top-level directory
    os.chdir(os.path.join(os.path.dirname(sys.argv[0]), '..'))
    unittest.main()
