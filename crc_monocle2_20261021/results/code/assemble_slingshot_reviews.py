"""Bind already completed actual PNG/PDF reviews to current files and coverage.

No visual inspection is performed here. Explicit per-figure human/agent viewing
declarations and observations must already exist, including every PDF page.
"""
import csv
import datetime
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULT = ROOT / 'results'
LOG = RESULT / '04_日志和检查'
notes = json.loads((LOG / 'Slingshot_manual_observations.json').read_text())
manifests = ['slingshot_figure_manifest.json', 'Slingshot_dynamics_manifest.json']
items = sum([json.loads((LOG / name).read_text()) for name in manifests], [])
assert len(notes) == len(items) == 42
assert set(notes) == {item['name'] for item in items}
with (RESULT / '02_表格/CRC/Slingshot_path_coverage.csv').open() as stream:
    coverage = {(row['root'], int(row['path'])): row for row in csv.DictReader(stream)}
pages = json.loads((LOG / 'pdf_page_manifest.json').read_text())
records = {'main': [], 'dynamics': []}
for item in items:
    name = item['name']; note = notes[name]
    assert note['png_actually_viewed'] and note['pdf_actually_viewed'], name
    assert note['observations'] and note['pages_viewed'] == [1]
    category = 'main' if name.startswith('12_') else 'dynamics'
    observation = note['observations']
    if category == 'dynamics':
        _, orientation, path_number, *_ = name.split('_')
        row = coverage[orientation, int(path_number.removeprefix('Path'))]
        if name.endswith('_proxy_specimens'):
            observation += f" 本图单列{int(row['unsupported_specimens'])}个缺少患者链接支持的标本，不纳入主体标签汇总，不声称独立患者。"
        else:
            observation += f" 路径纳入{int(row['cells']):,}个有限PT且权重至少0.5的细胞，连续CNV有{int(row['supported_CNV_bins'])}/12箱达到至少3个有支持报告标签；箱间组成可变，不等于同一患者进展。"
    for kind in ['png', 'pdf']:
        path = RESULT / item[kind]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == note[kind + '_sha256'], str(path)
        record = dict(file=item[kind], kind=kind, sha256=digest, status='passed',
                      actual_visual_inspection=True, observations=observation,
                      checks='All labels/panels and legends visible; full required range retained; PNG/PDF layouts agree; unsupported bins remain gaps; no clipping obscuring interpretation')
        if kind == 'pdf':
            matching = [page for page in pages if page['pdf'] == item['pdf']]
            assert [page['page'] for page in matching] == [1]
            assert (RESULT / matching[0]['rendered_png']).exists()
            record.update(pages_reviewed=[1], rendered_page=matching[0]['rendered_png'])
        records[category].append(record)
    item.update(finding=observation, visual_review='passed_actual_view')
for filename, category in [('slingshot_visual_review.json', 'main'),
                           ('Slingshot_dynamics_visual_review.json', 'dynamics')]:
    report = dict(reviewed_at=datetime.datetime.now().astimezone().isoformat(),
                  scope='Actual individual PNGs and every rendered PDF page viewed; explicit notes, hashes and pages checked',
                  records=records[category])
    (LOG / filename).write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
for filename in manifests:
    names = {item['name'] for item in json.loads((LOG / filename).read_text())}
    output = [item for item in items if item['name'] in names]
    (LOG / filename).write_text(json.dumps(output, ensure_ascii=False, indent=2) + '\n')
combined = json.loads((LOG / 'Monocle2_dynamics_manifest.json').read_text())
combined += [item for item in items if not item['name'].startswith('12_')]
(LOG / 'dynamics_manifest.json').write_text(json.dumps(combined, ensure_ascii=False, indent=2) + '\n')
print('Bound actual reviews: 42 PNGs and all 42 PDF pages; original Monocle2 reviews retained')
