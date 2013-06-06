import unittest
import sys
import shutil
import os
import subprocess
import modeller

class Tests(unittest.TestCase):
    def setUp(self):
        # Change into top-level directory and make directory for output
        os.chdir(os.path.join(os.path.dirname(sys.argv[0]), '..'))
        shutil.rmtree('output', ignore_errors=True)
        os.mkdir('output')

    def test_script1(self):
        """Test step 1 (build_profile)"""
        # Make sure the script runs without errors
        p = subprocess.check_call(['scripts/script1_build_profile.py'])
        # Make sure the profile contains the sequences we expect
        e = modeller.environ()
        a = modeller.alignment(e, file='output/build_profile.ali')
        self.assertEqual(sorted(s.code for s in a),
                         ['1a6dA', '1a6dB', '1dk7A', '1iokA', '1kidA',
                          '1kp8A', '1la1A', '1q3qA', '1sjpA', '1srvA',
                          '1we3A', '3kfbA', '3kfeA', 'P0A6F5'])
        os.unlink('output/build_profile.prf')
        os.unlink('output/build_profile.ali')

if __name__ == '__main__':
    unittest.main()
