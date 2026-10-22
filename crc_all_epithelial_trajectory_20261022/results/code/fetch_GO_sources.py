"""Freeze official public annotation files before any module enrichment."""
import datetime
import hashlib
import json
import re
import urllib.request

from pipeline_common import ROOT, LOG, save_json

sources = {
    'Homo_sapiens.gene_info.gz': 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz',
    'gene2go.gz': 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz',
    'go-basic.obo': 'https://current.geneontology.org/ontology/go-basic.obo',
    'NCBI_Gene_DATA_README.txt': 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/README',
}
folder = ROOT / 'data/GO_sources'
folder.mkdir(parents=True, exist_ok=True)
record = dict(status='downloading', files=[],
              sources_documentation=['https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/gene-processing/geneftp/',
                                     'https://www.geneontology.org/docs/download-ontology/',
                                     'https://www.geneontology.org/docs/go-annotations/'])
manifest_path = LOG / 'GO_source_downloads.json'
if manifest_path.exists():
    record = json.loads(manifest_path.read_text())
    record['status'] = 'resuming_verified_complete_files'
for name, url in sources.items():
    path = folder / name
    if path.exists():
        previous = [f for f in record['files'] if f['name'] == name]
        assert len(previous) == 1 and path.stat().st_size == previous[0]['bytes']
        check = hashlib.sha256()
        with path.open('rb') as f:
            for b in iter(lambda:f.read(8*1024**2),b''):check.update(b)
        assert check.hexdigest() == previous[0]['sha256']
        print('Reused verified complete annotation file',name,flush=True)
        continue
    if path.with_suffix(path.suffix + '.partial').exists():
        raise RuntimeError('Existing annotation source must be verified/versioned, never overwritten: ' + name)
    print('Downloading', name, flush=True)
    h = hashlib.sha256()
    temporary = path.with_suffix(path.suffix + '.partial')
    request_headers = {'User-Agent': 'Mozilla/5.0'}
    if name == 'go-basic.obo':
        request_headers['Range'] = 'bytes=0-'
    request = urllib.request.Request(url, headers=request_headers)
    with urllib.request.urlopen(request, timeout=60) as response, temporary.open('xb') as output:
        headers = {k: response.headers[k] for k in ['Last-Modified', 'ETag', 'Content-Length', 'Content-Type', 'Content-Range'] if k in response.headers}
        expected_length = int(headers['Content-Length']) if 'Content-Length' in headers else None
        if response.status == 206:
            match = re.fullmatch(r'bytes 0-(\d+)/(\d+)', headers.get('Content-Range',''))
            assert match and int(match[1])+1 == int(match[2]), 'Partial HTTP response does not contain the entire artifact'
            expected_length = int(match[2])
        for b in iter(lambda: response.read(1024**2), b''):
            output.write(b)
            h.update(b)
    assert expected_length is None or temporary.stat().st_size == expected_length
    temporary.replace(path)
    record['files'].append(dict(name=name, source_url=url, bytes=path.stat().st_size,
                               sha256=h.hexdigest(), selected_HTTP_headers=headers,
                               retrieved_at=datetime.datetime.now().astimezone().isoformat()))
    save_json(LOG / 'GO_source_downloads.json', record)
    print('Saved', name, path.stat().st_size, 'bytes', flush=True)
record['status'] = 'all_downloads_complete_pending_annotation_validation'
save_json(LOG / 'GO_source_downloads.json', record)
