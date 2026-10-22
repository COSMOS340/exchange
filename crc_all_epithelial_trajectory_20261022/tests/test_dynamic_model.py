"""Independent full dummy-design reference and aggregation stress fixtures."""
import sys
from pathlib import Path
import unittest

import numpy as np
import pandas as pd
from scipy.stats import f as f_distribution

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'results/code'))
from dynamic_model import spline_basis, bh_adjust, fit_dynamic_genes


def fixture():
    rng = np.random.default_rng(1201)
    unit = np.repeat(np.arange(40), 12)
    bins = np.tile(np.arange(12), 40)
    clusters = unit // 2
    obs = pd.DataFrame(dict(specimen_id=['S'+str(i) for i in unit],
                            cluster_unit=['P'+str(i) for i in clusters], bin=bins))
    x = spline_basis((bins+.5)/12)
    b = rng.normal(size=(5, 4))
    # Paired specimens share a bin-level nuisance component within each label.
    shared = rng.normal(scale=.4, size=(20, 12, 4))
    y = x @ b + rng.normal(scale=2, size=(40, 4))[unit] + shared[clusters, bins]
    return obs, y


class DynamicModelTests(unittest.TestCase):
    def test_coefficients_and_cluster_test_against_full_dummy_design(self):
        obs, y = fixture()
        actual = fit_dynamic_genes(obs, y, bins=12, chunk=2)
        units, _ = pd.factorize(obs.specimen_id, sort=False)
        clusters, unique = pd.factorize(obs.cluster_unit, sort=False)
        basis = spline_basis((obs.bin.to_numpy()+.5)/12)
        dummy = np.eye(40)[units]
        design = np.column_stack([dummy, basis])
        w = np.ones(len(obs))/12
        full_beta = np.linalg.lstsq(design*np.sqrt(w[:,None]), y*np.sqrt(w[:,None]), rcond=None)[0]
        np.testing.assert_allclose(actual['beta'], full_beta[-5:], atol=1e-11)
        bread = np.linalg.inv(design.T @ (w[:,None]*design))
        correction = 20/19 * (len(obs)-1)/(len(obs)-45)
        for gene in range(y.shape[1]):
            residual = y[:,gene]-design@full_beta[:,gene]
            meat = np.zeros((45,45))
            for group in range(len(unique)):
                keep = clusters==group
                score = design[keep].T @ (w[keep]*residual[keep])
                meat += np.outer(score,score)
            cov = (bread @ meat @ bread * correction)[-5:,-5:]
            beta = full_beta[-5:,gene]
            statistic = beta @ np.linalg.solve(cov,beta)/5
            self.assertAlmostEqual(actual['statistic'][gene],statistic,places=7)
            self.assertAlmostEqual(actual['p'][gene],f_distribution.sf(statistic,5,19),places=10)

    def test_unit_intercepts_and_gene_chunk_size_do_not_change_result(self):
        obs, y = fixture()
        unit, _ = pd.factorize(obs.specimen_id,sort=False)
        baseline = fit_dynamic_genes(obs,y,bins=12,chunk=1)
        changed = fit_dynamic_genes(obs,y+100*np.arange(40)[unit,None],bins=12,chunk=4)
        np.testing.assert_allclose(baseline['beta'],changed['beta'],atol=1e-10)
        np.testing.assert_allclose(baseline['p'],changed['p'],atol=1e-9)

    def test_correlated_specimens_are_not_independent_cells_or_specimens(self):
        obs, y = fixture()
        clustered = fit_dynamic_genes(obs,y,bins=12)
        obs['cluster_unit']=obs.specimen_id
        independent = fit_dynamic_genes(obs,y,bins=12)
        self.assertEqual(clustered['linkage_clusters'],20)
        self.assertEqual(independent['linkage_clusters'],40)
        self.assertTrue((clustered['statistic'] < independent['statistic']*.6).all())

    def test_constant_gene_and_duplicate_observation_handling(self):
        obs, y = fixture()
        out=fit_dynamic_genes(obs,np.ones((len(obs),1)),bins=12)
        self.assertTrue(np.isnan(out['p']).all())
        with self.assertRaises(ValueError):
            fit_dynamic_genes(pd.concat([obs,obs.iloc[:1]]),np.vstack([y,y[:1]]),bins=12)

    def test_BH_and_fixed_spline_domain(self):
        actual=bh_adjust([.01,.04,.03,.2,np.nan])
        np.testing.assert_allclose(actual[:4],[.04,4*.04/3,4*.04/3,.2])
        self.assertTrue(np.isnan(actual[-1]))
        self.assertEqual(spline_basis(np.array([0,.5,1])).shape,(3,5))
        with self.assertRaises(ValueError):
            spline_basis([1.01])


if __name__ == '__main__':
    unittest.main()
