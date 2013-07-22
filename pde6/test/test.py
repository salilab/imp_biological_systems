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

    def test_integrative_modeling_score(self):
        """Test the integrative modeling scoring"""
        # Make sure the script runs without errors
        p = subprocess.check_call(['./run_modeling.py', 'test'],
                                  cwd='integrative_modeling')

    def test_integrative_modeling(self):
        """Test the entire integrative modeling run"""
        import modeller
        # Compile the clustering program
        subprocess.check_call(['gfortran', 'cluster.f', 'u3best.f',
                               '-o', 'cluster.x'],
                              cwd='integrative_modeling/bin')

        # Run sampling
        subprocess.check_call(['./run_modeling.py'],
                              cwd='integrative_modeling')

        # Analysis
        subprocess.check_call(['bin/get_frames.sh'],
                              cwd='integrative_modeling')

        # Make sure that the three "known good" clusters are reproduced
        clusters = glob.glob('integrative_modeling/clustering/clus.*.pdb')
        clusters = [x for x in clusters if '-' not in x]
        exp_clusters = glob.glob('model_refinement/cluster*/model.pdb')

        env = modeller.environ()
        n_cluster = 0
        rms = []
        cluster_match = [0] * len(clusters)
        exp_cluster_match = [0] * len(exp_clusters)
        # Get a matrix of RMSD between all clusters and the expected clusters
        for ncluster, cluster in enumerate(clusters):
            per_cluster = []
            for nexp_cluster, exp_cluster in enumerate(exp_clusters):
                mc = modeller.model(env, file=cluster)
                s = modeller.selection(mc)
                a = modeller.alignment(env)
                me = modeller.model(env, file=exp_cluster)
                a.append_model(mc, align_codes='clus')
                a.append_model(me, align_codes='exp_clus')
                # We only care about the global (non-cutoff) RMSD, so use a
                # large cutoff so that refine_local doesn't increase the number
                # of equivalent positions at the expense of worsening the RMSD
                r = s.superpose(me, a, rms_cutoff=999.)
                if r.rms < 15.0:
                    cluster_match[ncluster] += 1
                    exp_cluster_match[nexp_cluster] += 1
                per_cluster.append(r.rms)
            rms.append(per_cluster)
        # Count the number of clusters which are close to an expected cluster
        ncluster_match = len(cluster_match) - cluster_match.count(0)
        # Count the number of expected clusters which are close to a cluster
        nexp_cluster_match = len(exp_cluster_match) - exp_cluster_match.count(0)
        # Make sure that each of the 3 expected clusters is close to one of
        # the clusters we produced (but not all the *same* cluster)
        self.assert_(ncluster_match >= 3 and nexp_cluster_match >= 3,
                     "Could not find any match between the %d clusters found "
                     "in this test and the 3 'known good' clusters (match "
                     "defined as all-atom RMSD less than 15.0A). "
                     "RMSD matrix: %s" % (len(clusters), str(rms)))

    def test_refinement(self):
        """Test the refinement script"""
        # Make sure the script runs without errors
        p = subprocess.check_call(['./model-single.py', '--run_quick_test'],
                                  cwd='model_refinement/cluster1')

if __name__ == '__main__':
    # Always run from top-level directory
    os.chdir(os.path.join(os.path.dirname(sys.argv[0]), '..'))
    unittest.main()
