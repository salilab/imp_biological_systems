#!/usr/bin/env python

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
        # Make directory for output
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

    def test_script5(self):
        """Test step 5 (template alignment)"""
        # Make sure the script runs without errors
        p = subprocess.check_call(['scripts/script5_template_alignment.py'])
        # Check output alignment
        e = modeller.environ()
        a = modeller.alignment(e, file='output/groel-1iokA.ali')
        self.assertEqual([x.code for x in a], ['1iok', 'P0A6F5'])
        os.unlink('output/groel-1iokA.ali')

    def test_script6(self):
        """Test step 6 (model building and assessment)"""
        # Get inputs (outputs from steps 3 and 5)
        shutil.copy('precalculate_results/stage3_density_segmentation/' \
                    'groel_subunit_11.mrc', 'output')
        shutil.copy('precalculate_results/stage5_template_alignment/' \
                    'groel-1iokA.ali', 'output')
        # Make sure the script runs without errors
        p = subprocess.check_call(['scripts/' \
                                   'script6_model_building_and_assessment.py'])
        # Check output models
        e = modeller.environ()
        for i in range(1, 11):
            base = 'output/P0A6F5.B9999%04d' % i
            pdb = base + '.pdb'
            trunc_pdb = base + '.truncated.pdb'
            trunc_fit_pdb = base + '.truncated.fitted.pdb'
            m = modeller.model(e, file=pdb)
            self.assertEqual(len(m.residues), 548)
            m = modeller.model(e, file=trunc_pdb)
            self.assertEqual(len(m.residues), 524)
            m = modeller.model(e, file=trunc_fit_pdb)
            self.assertEqual(len(m.residues), 524)
            os.unlink(pdb)
            os.unlink(trunc_pdb)
            os.unlink(trunc_fit_pdb)
        scores = 'output/model_building.scores.output'
        wc = len(open(scores).readlines())
        # Should be one line for each of 10 models, plus a header
        self.assertEqual(wc, 11)
        os.unlink(scores)

    def test_script7(self):
        """Test step 7 (pairwise RMSD)"""
        # Get inputs (outputs from steps 6)
        for i in range(1, 11):
            shutil.copy('precalculate_results/' \
                        'stage6_model_building_and_assessment/' \
                        'P0A6F5.B9999%04d.pdb' % i, 'output')
        # Make sure the script runs without errors
        p = subprocess.check_call(['scripts/' \
                                   'script7_pairwise_rmsd.py'])

    def test_script8(self):
        """Test step 8 (split density)"""
        # Make sure the script runs without errors
        p = subprocess.check_call(['scripts/script8_split_density.py'])
        os.unlink('output/groel-11.5A.top.mrc')
        os.unlink('output/groel-11.5A.bottom.mrc')

    def test_script9(self):
        """Test step 9 (multiple fitting)"""
        # Get inputs (outputs from step 8)
        for i in ('top', 'bottom'):
            shutil.copy('precalculate_results/stage8_split_density/' \
                        'groel-11.5A.%s.mrc' % i, 'output')
        # Make sure the script runs without errors
        p = subprocess.check_call(['scripts/' \
                                   'script9_symmetric_multiple_fitting.py'])
        e = modeller.environ()
        ref = modeller.model(e,
               file='precalculate_results/stage9_symmetric_multiple_fitting/' \
                    'model.top.0.pdb')
        sel = modeller.selection(ref).only_atom_types('CA')
        # At least one model in each ring should be close to the reference
        for side in ('top', 'bottom'):
            rms = []
            for i in range(6):
                fname = 'output/model.%s.%d.pdb' % (side, i)
                m =  modeller.model(e, file=fname)
                a = modeller.alignment(e)
                a.append_model(ref, align_codes='ref')
                a.append_model(m, align_codes='model')
                rms.append(sel.superpose(m, a).rms)
                os.unlink(fname)
            self.assert_(min(rms) < 10.0)
        os.unlink('output/intermediate_asmb_sols.out')
        for side in ('top', 'bottom'):
            os.unlink('output/multifit.%s.output' % side)
            os.unlink('output/multifit.%s.output.symm.ref' % side)
            os.unlink('output/multifit.%s.param' % side)

if __name__ == '__main__':
    # Always run from top-level directory
    os.chdir(os.path.join(os.path.dirname(sys.argv[0]), '..'))
    unittest.main()
