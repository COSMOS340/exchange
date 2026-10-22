"""Render every page of each existing figure PDF in a separate audit environment."""
import json
import hashlib
from pathlib import Path
import pymupdf as fitz

ROOT=Path(__file__).resolve().parents[2]
figures=ROOT/'results/01_图形'
destination=ROOT/'results/04_日志和检查/pdf_rendered'
manifest_file=ROOT/'results/04_日志和检查/pdf_page_manifest.json'
previous=json.loads(manifest_file.read_text()) if manifest_file.exists() else []
cache={(r['pdf'],r['page']):r for r in previous}
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
records=[]
new_pages=0
for pdf in sorted(figures.rglob('*.pdf')):
    relative=pdf.relative_to(figures)
    folder=destination/relative.parent
    folder.mkdir(parents=True,exist_ok=True)
    source_sha=digest(pdf)
    with fitz.open(pdf) as document:
        for number,page in enumerate(document,1):
            target=folder/f'{pdf.stem}_page{number:03d}.png'
            key=(str(pdf.relative_to(ROOT/'results')),number)
            old=cache.get(key)
            if old and old.get('source_pdf_sha256')==source_sha and target.exists() and old.get('rendered_png_sha256')==digest(target):
                records.append(old);continue
            page.get_pixmap(matrix=fitz.Matrix(150/72,150/72),alpha=False).save(target)
            new_pages+=1
            records.append(dict(pdf=str(pdf.relative_to(ROOT/'results')),page=number,
                rendered_png=str(target.relative_to(ROOT/'results')),
                source_pdf_sha256=source_sha,rendered_png_sha256=digest(target),
                total_pages=len(document),text_characters=len(page.get_text()),visual_review='pending'))
manifest_file.write_text(json.dumps(records,indent=2)+'\n')
print(f'Inventoried {len(records)} PDF pages; {new_pages} newly rendered. Actual visual review remains separate.')
