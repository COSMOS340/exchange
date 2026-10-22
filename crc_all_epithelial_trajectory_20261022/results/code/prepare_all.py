"""Prepare the full frozen matrix; no cell selection and no fitted models."""
import json
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from pipeline_common import (ROOT, BLOCK, CELLS, GENES, TAB, LOG, all_cell_mask,
                             map_sources, verify_data, digest, save_json, manifest, check_raw_counts)


def main():
    verify_data()
    work = ROOT / 'work'
    if work.exists():
        raise RuntimeError('New preparation requires an absent work directory')
    source = ROOT / 'data/Overall'
    original = pd.read_csv(source / 'metadata.csv.gz', low_memory=False)
    if len(original) != CELLS or not all_cell_mask(original).all():
        raise ValueError('Full-cell scope mismatch')
    scope = json.loads((ROOT / 'scope_summary.json').read_text())
    for field, expected in [('tissue_group', scope['tissue_group_2class'])]:
        if original[field].value_counts().to_dict() != expected:
            raise ValueError('Frozen counts differ: ' + field)
    audit_file = ROOT / 'prior_specimen_source_audit.csv'
    if digest(audit_file) != '71b626ffe5e42c051ae27d293ea87a0b394814182e353943c70b17a4c93542e8':
        raise ValueError('Frozen source audit SHA256 mismatch')
    m = map_sources(original, pd.read_csv(audit_file))
    if m.display_3class.value_counts().to_dict() != scope['tissue_source_3class']:
        raise ValueError('Three-class counts differ')
    units = dict(datasets=int(m.dataset.nunique()), specimen_labels=int(m.specimen_id.nunique()),
                 capture_units=int(m.capture_id.nunique()), subject_labels=int(m.subject_key.nunique()))
    if units != scope['frozen_units']:
        raise ValueError('Frozen statistical unit counts differ')
    for field, parents in [('capture_id', ['dataset', 'specimen_id', 'subject_key']),
                           ('specimen_id', ['dataset', 'subject_key'])]:
        if m.groupby(field)[parents].nunique().gt(1).any().any():
            raise ValueError('Unit hierarchy conflict: ' + field)
    genes = pd.read_csv(source / 'genes.csv')
    if len(genes) != GENES or not genes.gene.is_unique or genes.gene.isna().any():
        raise ValueError('Frozen gene identifiers differ')
    if genes.ordering.isna().any() or not genes.ordering.isin([True, False]).all():
        raise ValueError('Invalid ordering-gene flags')
    print('Exact source mapping passed for all cells; checking complete sparse matrix.', flush=True)
    with h5py.File(source / 'counts.h5') as h:
        shape = tuple(map(int, h['shape'][:]))
        if shape != (GENES, CELLS):
            raise ValueError('Frozen matrix dimensions differ')
        x = sp.csc_matrix((h['data'][:], h['indices'][:], h['indptr'][:]), shape=shape)
    detected = check_raw_counts(x, m, genes.gene)
    original_ordering = genes.ordering.to_numpy(dtype=bool)
    genes['source_ordering'] = original_ordering
    genes['source_detected_cells'] = genes.detected_cells
    if not np.array_equal(genes.detected_cells.to_numpy(), detected):
        raise ValueError('Source full-pool gene detection counts differ')
    genes['detected_cells'] = detected
    genes['ordering'] = original_ordering & (detected >= 10)
    if genes.ordering.sum() < 2:
        raise ValueError('Insufficient frozen ordering genes')
    nnz = int(x.nnz)
    del x
    target = work / BLOCK
    target.mkdir(parents=True)
    # Keep every original CSV field literally, including the CNV decimal text.
    # The numeric frame above validates values; this frame preserves delivery.
    literal = pd.read_csv(source / 'metadata.csv.gz', dtype=str, keep_default_na=False)
    if not literal.cell_id.equals(m.cell_id):
        raise ValueError('Literal metadata row order differs')
    for name in ['normal_origin', 'origin_evidence', 'display_3class']:
        literal[name] = m[name].to_numpy()
    literal.to_csv(target / 'metadata.csv.gz', index=False)
    del literal
    genes.to_csv(target / 'genes.csv', index=False)
    shutil.copy2(source / 'counts.h5', target / 'counts.h5')
    for name in ['Harmony', 'PCA']:
        z = np.load(source / (name + '.npy'), allow_pickle=False)
        if z.shape != (CELLS, 30) or not np.isfinite(z).all():
            raise ValueError('Invalid full-pool embedding: ' + name)
        pd.DataFrame(z).to_csv(target / (name + '.csv'), index=False, float_format='%.17g')
    cfg = json.loads((ROOT / 'config.json').read_text())
    cfg['blocks'] = [dict(block=BLOCK, identity='All frozen colorectal epithelial cells',
                          cells=CELLS, root='DEFERRED_TO_COMMON_GRAPH',
                          alternative='DEFERRED_TO_COMMON_GRAPH')]
    save_json(work / 'analysis_request.json', cfg)
    subprocess.run([sys.executable, str(ROOT / 'results/code/build_engine.py')], check=True)
    mapping = m[['cell_id', 'dataset', 'specimen_id', 'capture_id', 'subject_key',
                 'tissue_group', 'normal_origin', 'display_3class', 'origin_evidence']]
    mapping.to_csv(TAB / 'source_mapping_all_cells.csv.gz', index=False)
    for name in ['dataset', 'capture_id', 'specimen_id']:
        m.groupby([name, 'tissue_group', 'display_3class'], observed=True).size().rename('cells').reset_index().to_csv(
            TAB / ('source_counts_by_' + name + '.csv'), index=False)
    shutil.copy2(source / 'metadata.csv.gz', TAB / 'frozen_source_metadata.csv.gz')
    shutil.copy2(source / 'genes.csv', TAB / 'frozen_source_genes.csv')
    summary = dict(passed=True, input_cells=CELLS, fitted_cells_target=CELLS, sampled_cells=0,
                   genes=GENES, nnz=nnz, ordering_genes=int(genes.ordering.sum()),
                   tissue_counts=m.tissue_group.value_counts().to_dict(),
                   three_class_counts=m.display_3class.value_counts().to_dict(),
                   units=units, unmatched_source_cells=0, original_fields_unchanged=True,
                   UMI_and_cell_order_verified=True, all_cells_one_shared_tree=True,
                   roots='Selected after the common tree using endpoint/specimen support; labels unavailable to selector',
                   preparation_only=True)
    save_json(work / 'preparation_summary.json', summary)
    save_json(ROOT / 'WORK_MANIFEST.json', manifest(work))
    save_json(LOG / 'preparation_validation.json', summary)
    save_json(LOG / 'source_binding.json', dict(
        input_manifest_sha256=digest(ROOT / 'INPUT_MANIFEST.json'),
        source_audit_sha256=digest(audit_file),
        work_manifest_sha256=digest(ROOT / 'WORK_MANIFEST.json'),
        config_sha256=digest(ROOT / 'config.json'),
        analysis_plan_sha256=digest(ROOT / 'ANALYSIS_PLAN.json'),
        previous_CRC_result_commit='4d8a13008e57dc5a30a6fb5b5fc1a48410151cce'))
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
