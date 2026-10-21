"""Record the delivered source/environment and verify the original saved tree bytes."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
LOG = ROOT / 'results/04_日志和检查'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    fit = json.loads((LOG / 'Slingshot_local_process.json').read_text())
    assert fit['returncode'] == 0 and fit['status'] == 'fit_saved_pending_acceptance'
    assert (LOG / 'CRC_Slingshot_completed.json').exists()
    packages = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'], text=True)
    (LOG / 'python_freeze_final.txt').write_text(packages)
    expression = """pkgs=c('monocle','DDRTree','igraph','slingshot','princurve','Rcpp','Matrix','DelayedMatrixStats','sparseMatrixStats');
v=vapply(pkgs,function(p)as.character(packageVersion(p)),character(1));
stopifnot(v['monocle']=='2.38.0',v['DDRTree']=='0.1.6',v['igraph']=='2.2.2',v['slingshot']=='2.18.0',v['princurve']=='2.1.6');
write.table(data.frame(package=names(v),version=unname(v)),file=commandArgs(TRUE)[1],sep=',',row.names=FALSE);
sessionInfo()"""
    session = subprocess.check_output([str(ROOT / '.venv/bin/Rscript'), '-e', expression,
                                       str(LOG / 'R_actual_package_versions_final.csv')],
                                      cwd=ROOT, text=True, stderr=subprocess.STDOUT)
    (LOG / 'R_final_package_session.txt').write_text(session)
    old = json.loads((LOG / 'Monocle2_saved_object_hashes.json').read_text())
    baseline = {item['path']: item for item in old['files']}
    objects = []
    for path in sorted((ROOT / 'results/03_对象/CRC').glob('*.rds')):
        key = str(path.relative_to(ROOT))
        item = dict(path=key, bytes=path.stat().st_size, sha256=digest(path))
        if key in baseline:
            assert item == baseline[key], key
        objects.append(item)
        print('Hashed saved object:', path.name, flush=True)
    assert {item['path'] for item in objects}.issuperset(baseline)
    assert len(objects) == 7
    files = [path for path in (ROOT / 'results/code').rglob('*')
             if path.is_file() and '__pycache__' not in path.parts]
    files += [ROOT / name for name in ['work/fit_generated.R', 'WORK_MANIFEST.json',
                                       'INPUT_MANIFEST.json', 'config.json', 'Monocle2_PC.py']]
    sources = {str(path.relative_to(ROOT)): digest(path) for path in sorted(files)}
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    branch = subprocess.check_output(['git', 'branch', '--show-current'], cwd=ROOT, text=True).strip()
    changes = subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True)
    report = dict(passed=True, recorded_at=datetime.datetime.now().astimezone().isoformat(),
                  git_commit=commit, git_branch=branch, tracked_worktree_status=changes,
                  Python=sys.version, Python_executable=sys.executable, OS=platform.platform(),
                  logical_CPUs=os.cpu_count(), actual_R_versions='R_actual_package_versions_final.csv',
                  R_session='R_final_package_session.txt', Python_packages='python_freeze_final.txt',
                  conda_metadata_caveat='igraph metadata is 2.2.3, actual installed and used R package is 2.2.2 from official source',
                  original_Monocle2_object_bytes_unchanged=True,
                  saved_objects=objects, source_sha256=sources)
    (LOG / 'final_provenance.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    version = dict(input_commit='4bb17d4ffc23912ad1cd82ad87f4760d5e0f4955',
                   local_result_commit=commit, local_result_branch=branch,
                   Slingshot_fit_commit=json.loads((LOG / 'Slingshot_binding.json').read_text())['git_commit'],
                   full_source_binding='final_provenance.json',
                   remote_push_performed=False, result_release_published=False)
    (LOG / 'code_version.json').write_text(json.dumps(version, indent=2) + '\n')
    print('Final environment, source and seven real model objects recorded', flush=True)


if __name__ == '__main__':
    main()
