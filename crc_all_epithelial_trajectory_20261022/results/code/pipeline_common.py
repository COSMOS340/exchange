"""Shared paths, immutable input checks and auditable output helpers."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BLOCK = 'AllEpithelial'
CELLS = 223736
GENES = 15635
TAB = ROOT / 'results/02_表格' / BLOCK
OBJ = ROOT / 'results/03_对象' / BLOCK
LOG = ROOT / 'results/04_日志和检查'
ENV = Path('/home/anran/hehe/exchange/crc_monocle2_20261021/.venv')


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(b)
    return h.hexdigest()


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def manifest(folder):
    return [dict(path=p.relative_to(folder).as_posix(), bytes=p.stat().st_size,
                 sha256=digest(p)) for p in sorted(folder.rglob('*')) if p.is_file()]


def verify(folder, records):
    folder = Path(folder).resolve()
    for item in records:
        p = (folder / item['path']).resolve()
        if not p.is_relative_to(folder) or not p.is_file():
            raise ValueError('Unsafe or absent manifest file: ' + item['path'])
        if p.stat().st_size != item['bytes'] or digest(p) != item['sha256']:
            raise ValueError('Manifest mismatch: ' + item['path'])


def verify_data():
    verify(ROOT / 'data', json.loads((ROOT / 'INPUT_MANIFEST.json').read_text()))


def all_cell_mask(metadata):
    """Validate every frozen cell; never select a tissue subset."""
    import numpy as np
    required = ['cell_id', 'dataset', 'specimen_id', 'capture_id', 'sample_id',
                'subject_key', 'tissue_group', 'trajectory_eligible',
                'cnv_summary_eligible', 'identity_review', 'known_CD45_enriched',
                'CNV_relative_band', 'cnv_abs_log2', 'nCount_input', 'common_gene_UMI',
                'coarse_state', 'broad_identity', 'reported_patient_link_supported']
    if not set(required).issubset(metadata):
        raise ValueError('Missing frozen metadata fields')
    for name in ['cell_id', 'dataset', 'specimen_id', 'capture_id', 'sample_id',
                 'subject_key', 'coarse_state', 'broad_identity']:
        if metadata[name].isna().any():
            raise ValueError('Missing key: ' + name)
    if not metadata.cell_id.is_unique:
        raise ValueError('Duplicate cell_id')
    for name in ['trajectory_eligible', 'cnv_summary_eligible', 'identity_review',
                 'known_CD45_enriched', 'reported_patient_link_supported']:
        if metadata[name].isna().any() or not metadata[name].isin([True, False]).all():
            raise ValueError('Invalid boolean: ' + name)
    if not metadata.trajectory_eligible.all() or not metadata.cnv_summary_eligible.all():
        raise ValueError('Frozen input includes ineligible cells; do not delete them')
    if metadata.identity_review.any() or metadata.known_CD45_enriched.any():
        raise ValueError('Frozen identity flags differ; do not delete cells')
    if not metadata.tissue_group.isin(['Normal_colorectal', 'Primary_CRC']).all():
        raise ValueError('Unexpected exact tissue labels')
    if not metadata.CNV_relative_band.isin(['Low', 'Intermediate', 'High']).all():
        raise ValueError('Unexpected exact CNV labels')
    for name in ['cnv_abs_log2', 'nCount_input', 'common_gene_UMI']:
        if not np.isfinite(metadata[name]).all():
            raise ValueError('Nonfinite input: ' + name)
    if not metadata.nCount_input.gt(0).all() or not metadata.common_gene_UMI.gt(0).all():
        raise ValueError('Nonpositive UMI')
    if not metadata.sample_id.equals(metadata.capture_id):
        raise ValueError('sample_id and capture_id differ: audit before proceeding')
    return np.ones(len(metadata), dtype=bool)


def map_sources(metadata, audit):
    """Exact audited many-to-one join; row order and all original values remain."""
    import numpy as np
    keys = ['dataset', 'specimen_id', 'tissue_group']
    if audit.duplicated(keys).any() or audit[keys + ['normal_origin']].isna().any().any():
        raise ValueError('Source audit has duplicate or absent keys')
    original = metadata.copy()
    m = metadata.copy()
    m['_frozen_row'] = np.arange(len(m))
    m = m.merge(audit[keys + ['normal_origin', 'origin_evidence']], on=keys,
                how='left', validate='many_to_one', sort=False)
    if m.normal_origin.isna().any() or not np.array_equal(m._frozen_row, np.arange(len(m))):
        raise ValueError('Unmatched source or changed frozen row order')
    labels = {'Normal_unspecified_distance': 'Normal', 'Explicit_adjacent': 'Adjacent',
              'Not_normal': 'CRC'}
    if not m.normal_origin.isin(labels).all():
        raise ValueError('Unexpected exact normal_origin values')
    if not m.normal_origin.eq('Not_normal').equals(m.tissue_group.eq('Primary_CRC')):
        raise ValueError('Audited source conflicts with tissue_group')
    m['display_3class'] = m.normal_origin.map(labels)
    m = m.drop(columns='_frozen_row')
    if not m[original.columns].equals(original):
        raise ValueError('Source join changed an original field')
    return m


def fixed_bins(values, n=30):
    import numpy as np
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all() or values.min() < 0 or values.max() > 1:
        raise ValueError('Bins require finite values in [0,1]')
    return np.minimum(np.floor(n * values).astype(int), n - 1)


def equal_specimen_bins(frame, features, minimum_cells=10):
    """Cell -> eligible capture/bin -> equally weighted specimen/bin."""
    keys = ['dataset', 'subject_key', 'specimen_id', 'capture_id', 'bin']
    grouped = frame.groupby(keys, observed=True, dropna=False)
    captures = grouped[features].mean().join(grouped.size().rename('cells')).reset_index()
    captures['eligible'] = captures.cells.ge(minimum_cells)
    grouped = captures[captures.eligible].groupby(keys[:-2] + ['bin'], observed=True, dropna=False)
    specimens = grouped[features].mean().join(grouped.cells.sum().rename('cells')).join(
        grouped.size().rename('captures')).reset_index()
    return captures, specimens


def check_raw_counts(x, metadata, gene_names):
    import numpy as np
    if x.shape != (len(gene_names), len(metadata)):
        raise ValueError('Count dimensions differ from frozen ordered metadata')
    if not x.has_canonical_format or not np.issubdtype(x.data.dtype, np.integer) or (x.data < 0).any():
        raise ValueError('Expected canonical nonnegative integer counts')
    if not np.array_equal(np.asarray(x.sum(axis=0)).ravel(), metadata.common_gene_UMI.to_numpy()):
        raise ValueError('UMI totals differ from ordered metadata')
    return np.bincount(x.indices[x.data > 0], minlength=len(gene_names))


def validate_dynamic_source_tables(statistics, binned, smoothed, n_bins):
    import pandas as pd
    if statistics.gene.isna().any() or not statistics.gene.is_unique:
        raise ValueError('Dynamic statistics need unique complete genes')
    expected = pd.MultiIndex.from_product([statistics.gene, range(n_bins)], names=['gene', 'bin'])
    for name, table in [('binned', binned), ('smoothed', smoothed)]:
        if table[['gene', 'bin']].isna().any().any() or table.duplicated(['gene', 'bin']).any():
            raise ValueError('Duplicate or missing dynamic source keys: ' + name)
        observed = pd.MultiIndex.from_frame(table[['gene', 'bin']])
        if not observed.sort_values().equals(expected.sort_values()):
            raise ValueError('Dynamic source lacks an explicit gene/bin, including empty bins: ' + name)
    return dict(genes=len(statistics), bins=n_bins, rows_per_source=len(expected))
