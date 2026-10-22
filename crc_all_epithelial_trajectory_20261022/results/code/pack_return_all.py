"""Package accepted outputs, preserving objects, raw source tables and failures."""
import argparse,json,subprocess,zipfile
from pipeline_common import ROOT,LOG,digest,save_json


def included(path):
    return path.is_file() and '__pycache__' not in path.parts and path.suffix not in ['.pyc','.pyo'] and path.name!='.DS_Store'


def pack_return(label):
    if not label or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in label):
        raise ValueError('Use a new simple ASCII return label')
    for name in ['numeric_acceptance.json','visual_audit.json']:
        assert json.loads((LOG/name).read_text())['passed'],name
    assert all((ROOT/'results'/name).is_file() for name in ['RUN_SUMMARY.md','01_逐图阅读指南.md','02_方法与限制.md'])
    dirty=subprocess.check_output(['git','status','--porcelain','--untracked-files=all','--',str(ROOT)],cwd=ROOT,text=True)
    if dirty:raise RuntimeError('Commit reviewable code/docs before packaging: '+dirty)
    code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    target=ROOT/'return'/label;target.mkdir(parents=True,exist_ok=False)
    files=sorted(p for p in (ROOT/'results').rglob('*') if included(p))
    manifest=[dict(path=p.relative_to(ROOT/'results').as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in files]
    save_json(target/'RESULT_MANIFEST.json',manifest)
    archive=target/'results.zip'
    sidecars=['WORK_MANIFEST.json','INPUT_MANIFEST.json','config.json','ANALYSIS_PLAN.json','AllEpithelium_PC.py',
              'requirements.txt','TASK_SPEC.md','SOURCE_TISSUE_AUDIT.md','REFERENCE_FIGURE_REQUIREMENTS.md','RETURN_PROTOCOL.md','RUN_ON_PC.md','README.md','AGENTS.md']
    # Include every frozen source audit/manifest supplied in the task root.
    extras={ROOT/x for x in sidecars if (ROOT/x).is_file()}
    extras.update(p for p in ROOT.glob('*') if p.is_file() and p.suffix in ['.csv','.tsv','.json'])
    extras.update(p for p in (ROOT/'work').glob('*.json') if p.is_file())
    extras.add(ROOT/'work/fit_generated.R')
    sidecars=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(extras)]
    save_json(target/'SIDECAR_MANIFEST.json',sidecars)
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1,allowZip64=True) as zipped:
        for p in files:zipped.write(p,p.relative_to(ROOT).as_posix())
        for p in sorted(extras):zipped.write(p,p.relative_to(ROOT).as_posix())
        zipped.write(target/'RESULT_MANIFEST.json','RESULT_MANIFEST.json')
        zipped.write(target/'SIDECAR_MANIFEST.json','SIDECAR_MANIFEST.json')
    # Fail if any asynchronous process changed a result during packaging.
    for p,item in zip(files,manifest):assert p.stat().st_size==item['bytes'] and digest(p)==item['sha256'],item['path']
    cap=1536*1024**2;assets=[]
    if archive.stat().st_size<=cap:assets=[archive]
    else:
        with archive.open('rb') as source:
            number=1
            while True:
                block=source.read(8*1024**2)
                if not block:break
                part=target/f'results.zip.part{number:03d}'
                with part.open('xb') as output:
                    written=0
                    while block:
                        output.write(block);written+=len(block)
                        if written==cap:break
                        block=source.read(min(8*1024**2,cap-written))
                assets.append(part);number+=1
    info=dict(label=label,code_commit=code_commit,archive_bytes=archive.stat().st_size,archive_sha256=digest(archive),
        assets=[dict(name=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in assets],
        transfer_splitting_only=True,all_measured_cells=223736,original_input_archive_reuploaded=False,
        scientific_acceptance='PC numeric and actual visual checks passed; separate Mac receiving acceptance remains pending',
        exclusions=['Recomputable __pycache__/pyc/pyo and .DS_Store only; raw input and environment were never in results'])
    save_json(target/'RETURN_MANIFEST.json',info)
    print(json.dumps(info,indent=2),flush=True)
    return target


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--label',required=True);pack_return(p.parse_args().label)
