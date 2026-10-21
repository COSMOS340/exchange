"""Transcribe explicit actual PDF page reviews, never perform or infer a review.

Call only for the named figures whose rendered pages were actually viewed.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOG = ROOT / 'results/04_日志和检查'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('names', nargs='+')
args = parser.parse_args()
notes_path = LOG / 'Slingshot_manual_observations.json'
notes = json.loads(notes_path.read_text())
pages = json.loads((LOG / 'pdf_page_manifest.json').read_text())
for name in args.names:
    note = notes[name]
    assert note['png_actually_viewed'] and note['observations']
    folder = 'main' if name.startswith('12_') else 'dynamics'
    png = ROOT / 'results/01_图形' / folder / (name + '.png')
    assert hashlib.sha256(png.read_bytes()).hexdigest() == note['png_sha256']
    pdf = png.with_suffix('.pdf')
    matching = [page for page in pages if page['pdf'] == str(pdf.relative_to(ROOT / 'results'))]
    assert [page['page'] for page in matching] == [1]
    assert (ROOT / 'results' / matching[0]['rendered_png']).exists()
    note.update(pdf_actually_viewed=True, pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),
                pages_viewed=[1])
notes_path.write_text(json.dumps(notes, ensure_ascii=False, indent=2) + '\n')
print('Explicitly recorded viewed PDFs:', len(args.names),
      '; total', sum(note['pdf_actually_viewed'] for note in notes.values()))
