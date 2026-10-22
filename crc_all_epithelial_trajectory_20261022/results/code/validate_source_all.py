"""Independent streaming literal-field and source-rule audit of every cell."""
import argparse
import csv
import gzip
import itertools
from collections import Counter
from pathlib import Path

from pipeline_common import ROOT, CELLS, TAB, LOG, digest, save_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--previous-crc-dir', type=Path, required=True)
    args = parser.parse_args()
    with (ROOT / 'prior_specimen_source_audit.csv').open() as f:
        rows = list(csv.DictReader(f))
    audit = {(r['dataset'], r['specimen_id'], r['tissue_group']): r for r in rows}
    assert len(audit) == len(rows)
    counts, seen, crc_ids = Counter(), set(), []
    with gzip.open(ROOT / 'data/Overall/metadata.csv.gz', 'rt') as a, gzip.open(
            ROOT / 'work/AllEpithelial/metadata.csv.gz', 'rt') as b:
        source, derived = csv.DictReader(a), csv.DictReader(b)
        fields = source.fieldnames
        assert all(field in derived.fieldnames for field in fields)
        for original, work in itertools.zip_longest(source, derived):
            assert original is not None and work is not None
            assert all(original[k] == work[k] for k in fields), original['cell_id']
            cell = original['cell_id']
            assert cell not in seen
            seen.add(cell)
            key = tuple(original[k] for k in ['dataset', 'specimen_id', 'tissue_group'])
            expected = audit[key]['normal_origin']
            assert work['normal_origin'] == expected
            label = {'Normal_unspecified_distance': 'Normal', 'Explicit_adjacent': 'Adjacent',
                     'Not_normal': 'CRC'}[expected]
            assert work['display_3class'] == label
            assert (label == 'CRC') == (original['tissue_group'] == 'Primary_CRC')
            counts[label] += 1
            if label == 'CRC':
                crc_ids.append(cell)
    assert len(seen) == CELLS and counts == {'Normal': 64654, 'Adjacent': 8165, 'CRC': 150917}
    with gzip.open(args.previous_crc_dir / 'work/CRC/metadata.csv.gz', 'rt') as f:
        previous = [row['cell_id'] for row in csv.DictReader(f)]
    assert crc_ids == previous
    for name, target in [('metadata.csv.gz', 'frozen_source_metadata.csv.gz'),
                         ('genes.csv', 'frozen_source_genes.csv')]:
        assert digest(ROOT / 'data/Overall' / name) == digest(TAB / target)
    report = dict(passed=True, cells=len(seen), original_fields_compared=fields,
        every_original_field_literal_unchanged=True, CNV_decimal_text_and_bands_exact=True,
        source_mapping_independently_rebuilt_for_all_cells=True,
        source_counts=dict(counts), unmatched_cells=0,
        CRC_subset_exactly_matches_previous_ordered_150917_ids=True,
        original_metadata_copy_sha256=digest(TAB / 'frozen_source_metadata.csv.gz'),
        work_metadata_sha256=digest(ROOT / 'work/AllEpithelial/metadata.csv.gz'))
    save_json(LOG / 'frozen_source_to_work_comparison.json', report)
    print('All', len(seen), 'cells: every original metadata field literally unchanged, exact source mapping, and old CRC ID alignment passed.', flush=True)


if __name__ == '__main__':
    main()
