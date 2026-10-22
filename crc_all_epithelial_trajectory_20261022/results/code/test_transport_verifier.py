"""Generated archive fixtures: valid transport and three independent corruptions."""
import hashlib,json,subprocess,sys,tempfile,zipfile
from pathlib import Path
from pipeline_common import ROOT,LOG,save_json,digest


def fixture(folder,variant):
    folder.mkdir()
    payload={'results/图形源表.csv':'bin,value\n0,1\n1,2\n'.encode(),'results/code/demo.py':b'print("generated fixture only")\n','config.json':b'{"generated":true}\n'}
    sha=lambda value:hashlib.sha256(value).hexdigest()
    results=[dict(path=k.removeprefix('results/'),bytes=len(v),sha256=sha(v)) for k,v in payload.items() if k.startswith('results/')]
    sidecars=[dict(path='config.json',bytes=len(payload['config.json']),sha256=sha(payload['config.json']))]
    if variant=='wrong_sidecar_SHA':sidecars[0]['sha256']='0'*64
    save_json(folder/'RESULT_MANIFEST.json',results);save_json(folder/'SIDECAR_MANIFEST.json',sidecars)
    archive=folder/'results.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for k,v in payload.items():z.writestr(k,v)
        if variant=='unmanifested_member':z.writestr('not_manifested.txt','must be rejected')
        for p in ['RESULT_MANIFEST.json','SIDECAR_MANIFEST.json']:z.write(folder/p,p)
    data=archive.read_bytes();middle=len(data)//2
    parts=[]
    for number,value in enumerate([data[:middle],data[middle:]],1):
        path=folder/f'results.zip.part{number:03d}';path.write_bytes(value)
        parts.append(dict(name=path.name,bytes=len(value),sha256=sha(value)))
    save_json(folder/'RETURN_MANIFEST.json',dict(archive_bytes=len(data),archive_sha256=sha(data),assets=parts))
    if variant=='corrupt_part':
        p=folder/parts[0]['name'];v=bytearray(p.read_bytes());v[7]^=1;p.write_bytes(v)


def main():
    parent=Path(tempfile.mkdtemp(prefix='generated_transport_fixture_',dir=ROOT/'.runtime'))
    reports=[]
    for variant in ['valid','corrupt_part','wrong_sidecar_SHA','unmanifested_member']:
        folder=parent/variant;fixture(folder,variant)
        process=subprocess.run([sys.executable,str(ROOT/'results/code/verify_return_archive.py'),str(folder)],capture_output=True,text=True)
        assert (process.returncode==0)==(variant=='valid'),variant+' '+process.stdout+process.stderr
        reports.append(dict(generated_fixture=variant,returncode=process.returncode,expected_result_observed=True))
    save_json(LOG/'transport_verifier_software_tests.json',dict(passed=True,tests=reports,
        synthetic_archives_only=True,real_return_archive_verified=False,
        verifier_source_sha256=digest(ROOT/'results/code/verify_return_archive.py')))
    print('Generated valid archive accepted; corrupt part, wrong sidecar SHA and unmanifested member rejected.',flush=True)


if __name__=='__main__':main()
