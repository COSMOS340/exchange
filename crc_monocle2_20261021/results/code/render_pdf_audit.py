"""Render every page of each existing figure PDF in a separate audit environment."""
import json
from pathlib import Path
import pymupdf as fitz

ROOT=Path(__file__).resolve().parents[2]
figures=ROOT/'results/01_图形'
destination=ROOT/'results/04_日志和检查/pdf_rendered'
records=[]
for pdf in sorted(figures.rglob('*.pdf')):
    relative=pdf.relative_to(figures)
    folder=destination/relative.parent
    folder.mkdir(parents=True,exist_ok=True)
    with fitz.open(pdf) as document:
        for number,page in enumerate(document,1):
            target=folder/f'{pdf.stem}_page{number:03d}.png'
            page.get_pixmap(matrix=fitz.Matrix(150/72,150/72),alpha=False).save(target)
            records.append(dict(pdf=str(pdf.relative_to(ROOT/'results')),page=number,
                rendered_png=str(target.relative_to(ROOT/'results')),
                text_characters=len(page.get_text()),visual_review='pending'))
(ROOT/'results/04_日志和检查/pdf_page_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
print(f'Rendered {len(records)} PDF pages; human visual review still pending')
