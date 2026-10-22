import sys
from pathlib import Path
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'results/code'))
from cluster_profiles import adjusted_rand,standardize_profiles,exact_silhouettes,cluster_all_profiles


class ClusterTests(unittest.TestCase):
    def test_partition_metric_is_label_invariant(self):
        self.assertEqual(adjusted_rand([0,0,1,1],[8,8,3,3]),1)
        self.assertAlmostEqual(adjusted_rand([0,0,1,1],[0,1,0,1]),-.5)

    def test_gene_Z_is_invariant_to_gene_intercept_and_positive_scale(self):
        a=np.arange(40).reshape(5,8)
        np.testing.assert_allclose(standardize_profiles(a),standardize_profiles(a*3+20))
        self.assertTrue((standardize_profiles(np.ones((4,8)))==0).all())

    def test_silhouette_matches_independent_pair_loop_on_all_rows(self):
        x=np.random.default_rng(7).normal(size=(8,4))
        partitions={2:np.array([0]*4+[1]*4),3:np.array([0]*3+[1]*3+[2]*2)}
        actual=exact_silhouettes(x,partitions,chunk=3)
        for k,labels in partitions.items():
            expected=[]
            for i in range(len(x)):
                within=[np.linalg.norm(x[i]-x[j]) for j in range(len(x)) if labels[j]==labels[i] and i!=j]
                a=np.mean(within)
                b=min(np.mean([np.linalg.norm(x[i]-x[j]) for j in range(len(x)) if labels[j]==c]) for c in range(k) if c!=labels[i])
                expected.append((b-a)/max(a,b))
            np.testing.assert_allclose(actual[k],expected,atol=1e-12)

    def test_all_gene_clustering_recovers_separated_patterns(self):
        rng=np.random.default_rng(11)
        labels=np.repeat(np.arange(3),20)
        centers=np.eye(3)
        x=centers[labels]+rng.normal(scale=.01,size=(60,3))
        a=cluster_all_profiles(x,3,20261022)
        b=cluster_all_profiles(x,3,20261022)
        self.assertEqual(adjusted_rand(labels,a['labels']),1)
        np.testing.assert_array_equal(a['labels'],b['labels'])


if __name__=='__main__':unittest.main()
