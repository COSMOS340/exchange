"""Freeze official public annotation files before any module enrichment."""
import datetime
import hashlib
import urllib.request

from pipeline_common import ROOT, LOG, save_json

sources = {
    'Homo_sapiens.gene_info.gz': 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz',
    'gene2go.gz': 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz',
    'go-basic.obo': 'https://purl.obolibrary.org/obo/go/go-basic.obo',
    'NCBI_Gene_DATA_README.txt': 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/README',
}
folder = ROOT / 'data/GO_sources'
folder.mkdir(parents=True, exist_ok=True)
record = dict(status='downloading', files=[],
              sources_documentation=['https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/gene-processing/geneftp/',
                                     'https://www.geneontology.org/docs/download-ontology/',
                                     'https://www.geneontology.org/docs/go-annotations/'])
for name, url in sources.items():
    path = folder / name
    if path.exists() or path.with_suffix(path.suffix + '.partial').exists():
        raise RuntimeError('Existing annotation source must be verified/versioned, never overwritten: ' + name)
    print('Downloading', name, flush=True)
    h = hashlib.sha256()
    temporary = path.with_suffix(path.suffix + '.partial')
    with urllib.request.urlopen(url, timeout=60) as response, temporary.open('xb') as output:
        headers = {k: response.headers[k] for k in ['Last-Modified', 'ETag', 'Content-Length', 'Content-Type'] if k in response.headers}
        for b in iter(lambda: response.read(1024**2), b''):
            output.write(b)
            h.update(b)
    temporary.replace(path)
    record['files'].append(dict(name=name, source_url=url, bytes=path.stat().st_size,
                               sha256=h.hexdigest(), selected_HTTP_headers=headers,
                               retrieved_at=datetime.datetime.now().astimezone().isoformat()))
    save_json(LOG / 'GO_source_downloads.json', record)
    print('Saved', name, path.stat().st_size, 'bytes', flush=True)
record['status'] = 'all_downloads_complete_pending_annotation_validation'
save_json(LOG / 'GO_source_downloads.json', record)
