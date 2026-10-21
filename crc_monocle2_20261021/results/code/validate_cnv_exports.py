"""Compare frozen CNV decimal text and bands/indicators without float parsing."""
import csv
from decimal import Decimal
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
with gzip.open(ROOT/'work/CRC/metadata.csv.gz', 'rt') as stream:
    rows = list(csv.DictReader(stream))
original = {row['cell_id']: (row['cnv_abs_log2'], row['CNV_relative_band']) for row in rows}
assert len(rows) == len(original) == 150917
records = []
for name in ['function_and_CNV_values.csv.gz', 'Monocle2_primary_cells.csv.gz',
             'Monocle2_alternative_cells.csv.gz', 'Slingshot_primary_cells.csv.gz',
             'Slingshot_alternative_cells.csv.gz']:
    cells = []
    maximum = Decimal(0)
    differences = 0
    with gzip.open(ROOT/'results/02_表格/CRC'/name, 'rt') as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        has_cnv = 'cnv_abs_log2' in fields
        has_band = 'CNV_relative_band' in fields
        indicators = all('CNV_'+band+'_fraction' in fields for band in ['Low','Intermediate','High'])
        for row in reader:
            cell = row['cell_id']; cells.append(cell)
            value, band = original[cell]
            if has_cnv:
                difference = abs(Decimal(value)-Decimal(row['cnv_abs_log2']))
                maximum = max(maximum, difference)
                differences += difference != 0
            if has_band:
                assert row['CNV_relative_band'] == band
            if indicators:
                for category in ['Low','Intermediate','High']:
                    assert Decimal(row['CNV_'+category+'_fraction']) == int(category == band)
    assert len(cells) == len(set(cells)) == 150917 and set(cells) == set(original)
    assert maximum == 0
    records.append(dict(file=name, cells=len(cells), exact_ids=True,
        continuous_CNV_exported=has_cnv,
        maximum_absolute_decimal_difference=str(maximum) if has_cnv else None,
        differing_CNV_values=differences if has_cnv else None,
        band_field_verified=has_band, three_band_indicators_verified=indicators,
        metadata_join='Original frozen metadata used by cell_id; native coordinate tables need not duplicate CNV'))
report = dict(passed=True, scope='CNV decimal text and actual exported band fields/indicators',
              source_to_work_validation='frozen_source_to_work_comparison.json', files=records)
(ROOT/'results/04_日志和检查/CNV_export_precision_review.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
