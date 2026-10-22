"""Full-pool expression aggregation after the common tree has been validated.

Count I/O is streamed; neither cell fitting nor gene testing is split into models.
"""
import json

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from pipeline_common import ROOT, BLOCK, CELLS, TAB, OBJ, LOG, fixed_bins, save_json


def add_grouped(destination, codes, values):
    unique, local = np.unique(codes, return_inverse=True)
    incidence = sp.csr_matrix((np.ones(len(local)), (local, np.arange(len(local)))),
                              shape=(len(unique), len(local)))
    destination[unique] += (incidence @ values).toarray()


def main():
    if not json.loads((LOG / 'Monocle2_independent_validation.json').read_text())['passed']:
        raise RuntimeError('Saved full-cell common tree must be validated first')
    if (LOG / 'dynamic_expression_aggregation.json').exists():
        raise RuntimeError('Expression aggregation exists; validate/version rather than overwrite')
    m = pd.read_csv(ROOT / 'work' / BLOCK / 'metadata.csv.gz', low_memory=False)
    genes = pd.read_csv(ROOT / 'work' / BLOCK / 'genes.csv')
    assert len(m) == CELLS and m.cell_id.is_unique
    bins = json.loads((ROOT / 'ANALYSIS_PLAN.json').read_text())['dynamic_genes']['bins']
    folder = TAB / 'dynamic_genes'
    folder.mkdir(exist_ok=True)
    layouts = {}
    for root in ['primary', 'alternative']:
        cells = pd.read_csv(TAB / ('Monocle2_' + root + '_cells.csv.gz'))
        assert np.array_equal(cells.cell_id, m.cell_id)
        pt = cells.PT_1.to_numpy()
        scaled = (pt - pt.min()) / np.ptp(pt)
        cell_bins = fixed_bins(scaled, bins)
        codes, keys = pd.factorize(pd.MultiIndex.from_arrays([m.capture_id, cell_bins]), sort=False)
        first = np.full(len(keys), len(m), dtype=int)
        np.minimum.at(first, codes, np.arange(len(m)))
        group = m.iloc[first][['dataset', 'capture_id', 'specimen_id', 'subject_key',
                               'reported_patient_link_supported']].reset_index(drop=True)
        group['bin'] = cell_bins[first]
        group['cells'] = np.bincount(codes, minlength=len(keys))
        group['eligible'] = group.cells.ge(10)
        layouts[root] = dict(codes=codes, groups=group,
                             sums=np.zeros((len(group), len(genes)), dtype=float))
        assignment = m[['cell_id', 'dataset', 'capture_id', 'specimen_id', 'subject_key',
                        'reported_patient_link_supported']].copy()
        assignment['raw_pseudotime'] = pt
        assignment['scaled_pseudotime'] = scaled
        assignment['bin'] = cell_bins
        assignment.to_csv(folder / ('Monocle2_' + root + '_cell_bins.csv.gz'), index=False)
    spec_codes, specimens = pd.factorize(m.specimen_id, sort=False)
    detected_per_specimen = np.zeros((len(specimens), len(genes)), dtype=float)
    detected = np.zeros(len(genes), dtype=np.int64)
    with h5py.File(ROOT / 'work' / BLOCK / 'counts.h5') as h:
        ptr = h['indptr'][:]
        for lo in range(0, CELLS, 2000):
            hi = min(lo + 2000, CELLS)
            a, b = int(ptr[lo]), int(ptr[hi])
            x = sp.csr_matrix((h['data'][a:b].astype(float), h['indices'][a:b],
                               ptr[lo:hi+1] - a), shape=(hi-lo, len(genes)))
            assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(), m.common_gene_UMI.iloc[lo:hi])
            positive = x.copy()
            positive.data = (positive.data > 0).astype(float)
            positive.eliminate_zeros()
            detected += np.asarray(positive.sum(axis=0)).ravel().astype(np.int64)
            add_grouped(detected_per_specimen, spec_codes[lo:hi], positive)
            values = x.multiply((1e4/m.nCount_input.iloc[lo:hi].to_numpy())[:, None]).tocsr()
            values.data = np.log1p(values.data)
            for layout in layouts.values():
                add_grouped(layout['sums'], layout['codes'][lo:hi], values)
            if lo % 20000 == 0:
                print('Aggregated expression cells', hi, '/', CELLS, flush=True)
    assert np.array_equal(detected, genes.detected_cells)
    coverage = pd.DataFrame(dict(gene=genes.gene, detected_cells=detected,
        detected_cell_fraction=detected/CELLS,
        specimens_with_at_least_5_positive_cells=(detected_per_specimen>=5).sum(axis=0)))
    coverage['detection_eligible'] = coverage.detected_cell_fraction.ge(.01) & coverage.specimens_with_at_least_5_positive_cells.ge(10)
    coverage.to_csv(folder / 'dynamic_gene_detection_coverage.csv', index=False)
    results = {}
    for root, layout in layouts.items():
        groups = layout['groups']
        capture_means = layout['sums']
        capture_means /= groups.cells.to_numpy()[:, None]
        eligible = groups[groups.eligible].reset_index(drop=True)
        keys = ['dataset', 'specimen_id', 'subject_key', 'reported_patient_link_supported', 'bin']
        spec_group, levels = pd.factorize(pd.MultiIndex.from_frame(eligible[keys]), sort=False)
        spec_table = levels.to_frame(index=False)
        spec_table.columns = keys
        counts = np.bincount(spec_group, minlength=len(levels))
        spec_mean = np.zeros((len(levels), len(genes)))
        np.add.at(spec_mean, spec_group, capture_means[groups.eligible.to_numpy()])
        spec_mean /= counts[:, None]
        spec_table['captures'] = counts
        spec_cells = np.zeros(len(levels), dtype=np.int64)
        np.add.at(spec_cells, spec_group, eligible.cells.to_numpy())
        spec_table['cells'] = spec_cells
        spec_table['cluster_unit'] = np.where(spec_table.reported_patient_link_supported,
            'reported_link::' + spec_table.subject_key, 'unknown_links_pooled::' + spec_table.dataset)
        groups.to_csv(folder / ('Monocle2_' + root + '_capture_bins.csv.gz'), index=False)
        spec_table.to_csv(folder / ('Monocle2_' + root + '_specimen_bins.csv.gz'), index=False)
        destination = OBJ / ('gene_expression_bins_' + root + '.h5')
        with h5py.File(destination, 'x') as h:
            h.create_dataset('capture_bin_means', data=capture_means, compression='gzip', shuffle=True)
            h.create_dataset('specimen_bin_means', data=spec_mean, compression='gzip', shuffle=True)
            h.create_dataset('gene', data=genes.gene.to_numpy(dtype=object), dtype=h5py.string_dtype())
            h.create_dataset('capture_bin_cells', data=groups.cells.to_numpy())
        results[root] = dict(cells=CELLS, observed_capture_bins=len(groups),
                             eligible_capture_bins=int(groups.eligible.sum()),
                             observed_specimen_bins=len(spec_table), genes=len(genes),
                             specimen_clusters=int(spec_table.cluster_unit.nunique()))
        del layout['sums']
    save_json(LOG / 'dynamic_expression_aggregation.json', dict(
        completed=True, independently_validated=False, cells=CELLS, roots=results,
        expression='log1p(CP10K) using frozen nCount_input',
        units='Eligible capture-bin mean, then equal capture means in each specimen-bin',
        all_genes_detection_checked=True, count_streaming_is_not_model_splitting=True))
    print(json.dumps(results, indent=2), flush=True)


if __name__ == '__main__':
    main()
