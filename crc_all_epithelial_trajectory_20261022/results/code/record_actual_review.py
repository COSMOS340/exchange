"""Record an explicit observation only AFTER actually opening the named images.

This helper never inspects or approves an image by itself. The caller must
explicitly attest which PNG/PDF pages they opened and provide a concrete note.
"""
import argparse,json
from pipeline_common import ROOT,LOG,digest,save_json


def record(name,note,png_opened=False,pdf_pages=None):
    if len(note)<25:raise ValueError('Provide a concrete per-figure observation, not an automatic pass')
    files=list((ROOT/'results/01_图形').rglob(name+'.png'));assert len(files)==1
    png=files[0];pdf=png.with_suffix('.pdf');assert pdf.is_file()
    path=LOG/'actual_visual_observations.json'
    observations=json.loads(path.read_text()) if path.exists() else {}
    entry=observations.get(name,{})
    if entry.get('png_sha256')!=digest(png):entry={}
    if png_opened:entry.update(png_actually_viewed=True,png_sha256=digest(png))
    if pdf_pages:
        rendered=[r for r in json.loads((LOG/'pdf_page_manifest.json').read_text()) if r['pdf']==pdf.relative_to(ROOT/'results').as_posix()]
        assert rendered and all(r['source_pdf_sha256']==digest(pdf) for r in rendered)
        expected=set(range(1,rendered[0]['total_pages']+1));viewed=set(pdf_pages)
        if entry.get('pdf_sha256')==digest(pdf):viewed.update(entry.get('pdf_pages_viewed',[]))
        assert viewed<=expected
        entry.update(pdf_actually_viewed=viewed==expected,pdf_pages_viewed=sorted(viewed),pdf_sha256=digest(pdf))
    entry['observation']=note
    observations[name]=entry;save_json(path,observations)
    print('Explicit opened-image observation recorded:',name,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--note',required=True)
    p.add_argument('--png-opened',action='store_true');p.add_argument('--pdf-pages',nargs='+',type=int)
    a=p.parse_args();record(a.name,a.note,a.png_opened,a.pdf_pages)
