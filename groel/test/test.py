import unittest
import sys
import shutil
import os
import glob
import subprocess
import modeller

class Tests(unittest.TestCase):
    templates = ['1a6dA', '1a6dB', '1dk7A', '1iokA', '1kidA', '1kp8A',
                 '1la1A', '1q3qA', '1sjpA', '1srvA', '1we3A', '3kfbA', '3kfeA']

    def setUp(self):
        # Change into top-level directory and make directory for output
        os.chdir(os.path.join(os.path.dirname(sys.argv[0]), '..'))
        shutil.rmtree('output', ignore_errors=True)
        os.mkdir('output')

    def test_script1(self):
        """Test step 1 (build profile)"""
        # Make sure the script runs without errors
        p = subprocess.check_call(['scripts/script1_build_profile.py'])
        # Make sure the profile contains the sequences we expect
        e = modeller.environ()
        a = modeller.alignment(e, file='output/build_profile.ali')
        self.assertEqual(sorted(s.code for s in a),
                         sorted(self.templates) + ['P0A6F5'])
        os.unlink('output/build_profile.prf')
        os.unlink('output/build_profile.ali')

    def test_script2(self):
        """Test step 2 (compare templates)"""
        # Get inputs (outputs from step 1)
        shutil.copy('precalculate_results/stage1_template_identification/' \
                    'build_profile.prf', 'output')
        # Make sure the script runs without errors
        p = subprocess.check_call(['scripts/script2_compare_templates.py'])
        os.unlink('output/family.mat')

    def test_script3(self):
        """Test step 3 (density segmentation)"""
        # Make sure the script runs without errors
        p = subprocess.check_call(['scripts/script3_density_segmentation.py'])
        # Should have produced a PDB with coordinates of all 14 subunit centers
        e = modeller.environ()
        m = modeller.model(e, file='output/groel_segments_center.pdb')
        self.assertEqual(len(m.atoms), 14)
        self.assertEqual(len(m.residues), 14)
        # load_configuration file should load all 14 subunits, and set level
        wc = len(open('output/load_configuration.cmd').readlines())
        self.assertEqual(wc, 15)
        os.unlink('output/load_configuration.cmd')
        os.unlink('output/groel_segments_center.pdb')
        for i in range(14):
            os.unlink('output/groel_subunit_%d.mrc' % i)

    def test_script4(self):
        """Test step 4 (score templates by cc)"""
        # Get inputs (outputs from step 1 and 3)
        shutil.copy('precalculate_results/stage1_template_identification/' \
                    'build_profile.prf', 'output')
        shutil.copy('precalculate_results/stage3_density_segmentation/' \
                    'groel_subunit_11.mrc', 'output')
        # Make sure the script runs without errors
        p = subprocess.check_call(['scripts/script4_score_templates_by_cc.py'])
        # Check output PDBs
        e = modeller.environ()
        long_templates = self.templates[:]
        # Remove templates that weren't long enough
        for t in ['1dk7A', '1kidA', '1la1A', '1srvA']:
            long_templates.remove(t)
        for code in long_templates:
            fname = 'output/%s_fitted.pdb' % code
            m = modeller.model(e, file=fname)
            self.assert_(len(m.atoms) > 10)
            os.unlink(fname)
        os.unlink('output/score_templates_by_cc.log')

if __name__ == '__main__':
    unittest.main()
