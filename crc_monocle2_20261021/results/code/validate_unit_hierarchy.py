"""Check the frozen library/specimen/reported-label grain without changing it."""
import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
m = pd.read_csv(ROOT/'work/CRC/metadata.csv.gz', low_memory=False)
keys = ['cell_id','dataset','sample_id','specimen_id','subject_key','reported_patient_link_supported']
assert len(m) == 150917 and m.cell_id.is_unique and not m[keys].isna().any().any()
checks = []
for unit, fields in [
    (['dataset','sample_id'], ['specimen_id','subject_key','reported_patient_link_supported']),
    (['dataset','specimen_id'], ['subject_key','reported_patient_link_supported']),
    (['dataset','subject_key'], ['reported_patient_link_supported']),
]:
    counts = m.groupby(unit, observed=True)[fields].nunique(dropna=False)
    bad = counts.gt(1).any(axis=1)
    checks.append(dict(unit=unit, units=len(counts), mapped_fields=fields,
                       conflicting_units=int(bad.sum())))
    assert not bad.any(), counts[bad].to_dict()
cohorts = m.groupby('dataset', observed=True).agg(cells=('cell_id','size'),
    specimens=('specimen_id','nunique'), libraries=('sample_id','nunique'),
    reported_labels=('subject_key','nunique'),
    supported_link_cells=('reported_patient_link_supported','sum'))
cohorts.to_csv(ROOT/'results/02_表格/CRC/frozen_unit_hierarchy_by_cohort.csv')
report = dict(passed=True, cells=len(m), cohorts=m.dataset.nunique(),
    specimens=len(m[['dataset','specimen_id']].drop_duplicates()),
    libraries=len(m[['dataset','sample_id']].drop_duplicates()),
    reported_labels=len(m[['dataset','subject_key']].drop_duplicates()),
    supported_reported_labels=len(m.loc[m.reported_patient_link_supported,['dataset','subject_key']].drop_duplicates()),
    independent_patient_identities_verified=False,
    interpretation='Consistency of frozen identifiers only; no proof of independent patient identities',
    hierarchical_checks=checks)
(ROOT/'results/04_日志和检查/frozen_unit_hierarchy_validation.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
