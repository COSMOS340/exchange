"""Run one bound full-cell fit with a lock, actual PID and peak-memory log."""
import datetime
import fcntl
import json
import os
from pathlib import Path
import resource
import subprocess
import time

from pipeline_common import ROOT, BLOCK, CELLS, OBJ, TAB, LOG, ENV, digest, save_json, verify, verify_data


def now():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()


def main(method):
    for p in [LOG, TAB, OBJ, ROOT / '.runtime/tmp', ROOT / '.runtime/cache']:
        p.mkdir(parents=True, exist_ok=True)
    with (LOG / 'fit.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        record_file = LOG / (method + '_local_process.json')
        if record_file.exists():
            raise RuntimeError('Inspect and preserve the prior version before retrying')
        verify_data()
        verify(ROOT / 'work', json.loads((ROOT / 'WORK_MANIFEST.json').read_text()))
        if not json.loads((LOG / 'preparation_validation.json').read_text())['passed']:
            raise RuntimeError('Full-cell preparation is not validated')
        for name in ['graph_compatibility_validation.json', 'root_fallback_synthetic_validation.json',
                     'endpoint_root_synthetic_validation.json', 'python_unit_test_result.json']:
            if not json.loads((LOG / name).read_text())['passed']:
                raise RuntimeError('Software validation failed: ' + name)
        if method == 'Slingshot':
            if not json.loads((LOG / 'Monocle2_independent_validation.json').read_text())['passed']:
                raise RuntimeError('Validate the saved common tree before Slingshot')
            tested = json.loads((LOG / 'exact_projection_software_binding.json').read_text())
            for name, expected in tested['tested_source_sha256'].items():
                if digest(ROOT / name) != expected:
                    raise RuntimeError('Projection source changed since testing: ' + name)
            entry = ROOT / 'results/code/fit_slingshot_exact_entry.R'
        else:
            entry = ROOT / 'work/fit_generated.R'
        completed = LOG / (BLOCK + '_' + method + '_completed.json')
        if completed.exists() or any(OBJ.glob(method + '*.rds')):
            raise RuntimeError('Existing formal objects must be preserved and reviewed')
        env = os.environ.copy()
        env.pop('R_MAX_VSIZE', None)
        env['PATH'] = str(ENV / 'bin') + os.pathsep + env.get('PATH', '')
        env.update(TMPDIR=str(ROOT / '.runtime/tmp'), XDG_CACHE_HOME=str(ROOT / '.runtime/cache'),
                   MPLCONFIGDIR=str(ROOT / '.runtime/matplotlib'), TZ='Asia/Shanghai')
        for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']:
            env[key] = '2'
        preflight = """stopifnot(getRversion()>='4.5.0');
pkgs=c('monocle','DDRTree','igraph','Rcpp','Matrix','data.table','jsonlite','rhdf5','dplyr','rlang','R.utils','slingshot','princurve','DelayedMatrixStats','sparseMatrixStats');
stopifnot(all(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)));
stopifnot(as.character(packageVersion('monocle'))=='2.38.0',as.character(packageVersion('DDRTree'))=='0.1.6',as.character(packageVersion('igraph'))=='2.2.2',as.character(packageVersion('slingshot'))=='2.18.0',as.character(packageVersion('princurve'))=='2.1.6');
Rcpp::cppFunction('int compiler_check(){return 1;}');stopifnot(compiler_check()==1);sessionInfo()"""
        pre = subprocess.run([str(ENV / 'bin/Rscript'), '-e', preflight], env=env, text=True, capture_output=True)
        (LOG / (method + '_preflight.txt')).write_text(pre.stdout + pre.stderr)
        if pre.returncode:
            raise RuntimeError('Actual environment preflight failed; inspect retained log')
        code_sources = [entry, ROOT / 'results/code/run_fit_all.py', ROOT / 'results/code/select_endpoint_roots.R']
        if method == 'Slingshot':
            code_sources += [ROOT / 'results/code/exact_projection_acceleration.R', ROOT / 'results/code/exact_projection_bvh.cpp']
        binding = dict(method=method, cells=CELLS, code_commit=subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            work_manifest_sha256=digest(ROOT / 'WORK_MANIFEST.json'),
            analysis_plan_sha256=digest(ROOT / 'ANALYSIS_PLAN.json'),
            source_sha256={p.relative_to(ROOT).as_posix(): digest(p) for p in code_sources},
            Rscript=str(ENV / 'bin/Rscript'), threads=2, original_environment_not_modified=True)
        save_json(LOG / (method + '_binding.json'), binding)
        command = [str(ENV / 'bin/Rscript'), str(entry), str(ROOT / 'work' / BLOCK),
                   str(OBJ), str(TAB), str(LOG), str(ROOT / 'work/analysis_request.json'), method]
        timed = ['/usr/bin/time', '-v', '-o', str(LOG / (method + '_resource_usage.txt')), *command]
        record = dict(method=method, command=command, cwd=str(ROOT), started_at=now(),
                      status='running', cells=CELLS, returncode=None,
                      formal_fit_started=True, analysis_accepted=False)
        started = time.monotonic()
        with (LOG / (method + '_console.log')).open('x') as console:
            child = subprocess.Popen(timed, cwd=ROOT, env=env, stdout=console,
                                     stderr=subprocess.STDOUT, start_new_session=True)
            record['supervisor_child_pid'] = child.pid
            save_json(record_file, record)
            save_json(LOG / 'progress.json', record)
            print(now(), method, 'started; all', CELLS, 'cells fitted jointly', flush=True)
            with (ROOT / 'results/RUN_SUMMARY.md').open('a') as f:
                f.write('\n' + now() + '：启动全量' + method + '，所有223,736细胞共同拟合；PID/命令/资源/退出状态记录于日志。保存完成不代表验收。\n')
            code = child.wait()
        success = code == 0 and completed.exists()
        record.update(status='fit_saved_pending_acceptance' if success else 'fit_failed',
                      returncode=code, finished_at=now(), elapsed_seconds=time.monotonic() - started,
                      children_peak_RSS_kB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
        save_json(record_file, record)
        save_json(LOG / (method + '_exit.json'), dict(returncode=code, completed_marker_exists=completed.exists()))
        save_json(LOG / 'progress.json', record)
        with (ROOT / 'results/RUN_SUMMARY.md').open('a') as f:
            f.write('\n' + now() + '：' + method + '退出码' + str(code) + '，耗时' + str(round(record['elapsed_seconds'], 1)) + '秒；状态' + record['status'] + '。\n')
        print(json.dumps(record, ensure_ascii=False, indent=2), flush=True)
        return 0 if success else 1
