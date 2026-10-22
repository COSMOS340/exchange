"""Scope and scientific aggregation regression tests; no biological models."""
import sys
from pathlib import Path
import unittest

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'results/code'))
from pipeline_common import (all_cell_mask, map_sources, fixed_bins, equal_specimen_bins,
                             check_raw_counts, validate_dynamic_source_tables)


def metadata(n=20):
    ids = np.arange(n)
    return pd.DataFrame(dict(cell_id=['v' + str(i) for i in ids], dataset='D',
        specimen_id=np.where(ids % 2 == 0, 'S0', 'S1'), capture_id='L', sample_id='L',
        subject_key='P', tissue_group=np.where(ids % 2 == 0, 'Normal_colorectal', 'Primary_CRC'),
        trajectory_eligible=True, cnv_summary_eligible=True, identity_review=False,
        known_CD45_enriched=False, CNV_relative_band='Low', cnv_abs_log2=.01,
        nCount_input=20, common_gene_UMI=10, coarse_state='C0', broad_identity='Epithelial',
        reported_patient_link_supported=True))


def audit():
    return pd.DataFrame(dict(dataset=['D', 'D'], specimen_id=['S0', 'S1'],
        tissue_group=['Normal_colorectal', 'Primary_CRC'],
        normal_origin=['Normal_unspecified_distance', 'Not_normal'], origin_evidence=['frozen', 'frozen']))


class FullCellTests(unittest.TestCase):
    def test_all_223736_cells_retained_including_normal(self):
        m = metadata(223736)
        mask = all_cell_mask(m)
        self.assertTrue(mask.all())
        self.assertEqual(mask.sum(), 223736)
        self.assertEqual(m.loc[mask, 'tissue_group'].value_counts()['Normal_colorectal'], 111868)

    def test_exact_three_class_join_preserves_order_and_values(self):
        m = metadata()
        a = audit()
        a.loc[0, 'normal_origin'] = 'Explicit_adjacent'
        result = map_sources(m, a.iloc[::-1])
        self.assertTrue(result[m.columns].equals(m))
        self.assertEqual(result.display_3class.tolist(), ['Adjacent', 'CRC'] * 10)

    def test_source_key_mismatch_stops_instead_of_name_guessing(self):
        a = audit(); a.loc[0, 'specimen_id'] = 'S0-N'
        with self.assertRaises(ValueError):
            map_sources(metadata(), a)
        with self.assertRaises(ValueError):
            map_sources(metadata(), pd.concat([audit(), audit().iloc[:1]]))

    def test_ineligible_duplicate_and_nonfinite_are_not_deleted(self):
        for field, value in [('trajectory_eligible', False), ('cnv_abs_log2', np.nan),
                             ('CNV_relative_band', 'low'), ('cell_id', 'v1')]:
            m = metadata(); m.loc[0, field] = value
            with self.assertRaises(ValueError):
                all_cell_mask(m)

    def test_source_tissue_conflict_rejected(self):
        a = audit(); a.loc[0, 'normal_origin'] = 'Not_normal'
        with self.assertRaises(ValueError):
            map_sources(metadata(), a)

    def test_matrix_order_and_integer_counts(self):
        x = sp.csc_matrix(np.array([[1, 0, 4], [0, 2, 1]], dtype=np.int32))
        m = pd.DataFrame({'common_gene_UMI': [1, 2, 5]})
        np.testing.assert_array_equal(check_raw_counts(x, m, ['A', 'B']), [2, 2])
        with self.assertRaises(ValueError):
            check_raw_counts(x[:, [2, 1, 0]], m, ['A', 'B'])
        with self.assertRaises(ValueError):
            check_raw_counts(x.astype(float), m, ['A', 'B'])

    def test_fixed_bin_endpoints_and_invalid_values(self):
        np.testing.assert_array_equal(fixed_bins([0, 1/30, .5, 1]), [0, 1, 15, 29])
        for bad in [[np.nan], [-.001], [1.001]]:
            with self.assertRaises(ValueError):
                fixed_bins(bad)

    def test_equal_capture_then_specimen_not_cell_weighting(self):
        frame = pd.DataFrame(dict(dataset='D', subject_key='P', specimen_id=['S1'] * 100 + ['S2'] * 10,
            capture_id=['A'] * 90 + ['B'] * 10 + ['C'] * 10, bin=0,
            expression=[10.] * 90 + [0.] * 10 + [1.] * 10))
        captures, specimens = equal_specimen_bins(frame, ['expression'])
        self.assertEqual(captures.cells.tolist(), [90, 10, 10])
        self.assertEqual(specimens.expression.tolist(), [5., 1.])
        self.assertEqual(specimens.expression.mean(), 3.)
        self.assertNotEqual(specimens.expression.mean(), frame.expression.mean())

    def test_empty_bins_are_explicit_and_duplicate_sources_rejected(self):
        statistics = pd.DataFrame({'gene': ['A', 'B']})
        table = pd.DataFrame({'gene': ['A', 'A', 'B', 'B'], 'bin': [0, 1, 0, 1],
                              'mean': [0., np.nan, 1., 2.]})
        self.assertEqual(validate_dynamic_source_tables(statistics, table, table, 2)['rows_per_source'], 4)
        with self.assertRaises(ValueError):
            validate_dynamic_source_tables(statistics, table.iloc[:-1], table, 2)
        with self.assertRaises(ValueError):
            validate_dynamic_source_tables(statistics, pd.concat([table, table.iloc[:1]]), table, 2)


if __name__ == '__main__':
    unittest.main()
