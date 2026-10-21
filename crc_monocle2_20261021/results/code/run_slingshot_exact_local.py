"""Versioned full-cell Slingshot run using the verified exact segment search.

Checks original input/work manifests and the tested runtime source binding.
Never modifies the frozen engine or re-fits Monocle2.
"""
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
LOG=ROOT/'results/04_日志和检查'
sys.path.insert(0,str(ROOT))
import Monocle2_PC as original


def now():
    return datetime.datetime.now().astimezone().isoformat()


def save(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
    temporary.replace(path)


def main():
    with (LOG/'fit.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        record_file=LOG/'Slingshot_local_process.json'
        if record_file.exists() or (LOG/'CRC_Slingshot_completed.json').exists():
            raise RuntimeError('Inspect and archive prior Slingshot attempt before a new version')
        if any((ROOT/'results/03_对象/CRC').glob('Slingshot*.rds')):
            raise RuntimeError('Existing Slingshot objects require preserved versioned handling')
        report=json.loads((LOG/'Monocle2_independent_validation.json').read_text())
        assert report['passed']
        software=json.loads((LOG/'exact_projection_synthetic_validation.json').read_text())
        assert software['passed'] and not software['formal_biological_acceleration_applied']
        tested=json.loads((LOG/'exact_projection_software_binding.json').read_text())
        for name,digest in tested['tested_source_sha256'].items():
            assert original.digest(ROOT/name)==digest, name
        original.verify_data()
        original.verify(ROOT/'work',json.loads((ROOT/'WORK_MANIFEST.json').read_text()))
        assert original.digest(ROOT/'work/fit_generated.R')=='91079e6bb6555671c35f3f114fdc64797e39d3b460c7f60fefebd4691ce85161'
        rscript=ROOT/'.venv/bin/Rscript'
        env=os.environ.copy();env.pop('R_MAX_VSIZE',None)
        for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
            env[key]='2'  # Same formal fit environment as original Monocle2_PC.fit.
        preflight="""stopifnot(getRversion()>='4.5.0');
pkgs=c('monocle','DDRTree','igraph','Rcpp','Matrix','data.table','jsonlite','rhdf5','dplyr','rlang','R.utils','slingshot','princurve','DelayedMatrixStats','sparseMatrixStats');
stopifnot(all(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)));
stopifnot(as.character(packageVersion('monocle'))=='2.38.0',as.character(packageVersion('DDRTree'))=='0.1.6',as.character(packageVersion('igraph'))=='2.2.2',as.character(packageVersion('slingshot'))=='2.18.0',as.character(packageVersion('princurve'))=='2.1.6');
Rcpp::cppFunction('int compiler_check(){return 1;}');stopifnot(compiler_check()==1);sessionInfo()"""
        checked=subprocess.run([str(rscript),'-e',preflight],env=env,capture_output=True,text=True)
        (LOG/'Slingshot_preflight.txt').write_text(checked.stdout+checked.stderr)
        if checked.returncode:raise RuntimeError('Environment preflight failed; inspect retained log')
        source_names=['work/fit_generated.R','results/code/fit_slingshot_exact_entry.R',
                      'results/code/exact_projection_acceleration.R','results/code/exact_projection_bvh.cpp',
                      'results/code/run_slingshot_exact_local.py']
        binding=dict(method='Slingshot',version='Slingshot_exact_search_v3',
            work_manifest_sha256=original.digest(ROOT/'WORK_MANIFEST.json'),
            Rscript=str(rscript),source_sha256={x:original.digest(ROOT/x) for x in source_names},
            preflight_sha256=hashlib.sha256(checked.stdout.encode()).hexdigest(),
            software_validation_sha256=original.digest(LOG/'exact_projection_synthetic_validation.json'),
            original_cells_and_curve_vertices_retained=True,fitting_parameters_changed=False,
            git_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
        save(LOG/'Slingshot_binding.json',binding)
        command=[str(rscript),str(ROOT/'results/code/fit_slingshot_exact_entry.R'),
                 str(ROOT/'work/CRC'),str(ROOT/'results/03_对象/CRC'),str(ROOT/'results/02_表格/CRC'),
                 str(LOG),str(ROOT/'work/analysis_request.json'),'Slingshot']
        record=dict(method='Slingshot',version='Slingshot_exact_search_v3',command=command,
                    cwd=str(ROOT),started_at=now(),status='running',returncode=None,
                    formal_fit_started=True,analysis_accepted=False)
        started=time.monotonic()
        with (LOG/'Slingshot_console.log').open('x') as console:
            child=subprocess.Popen(command,cwd=ROOT,env=env,stdout=console,stderr=subprocess.STDOUT,
                                   start_new_session=True)
            record['R_child_pid']=child.pid;save(record_file,record);save(LOG/'progress.json',record)
            print(now(),'Slingshot exact search v3 started; all 150917 cells and native parameters retained',flush=True)
            with (ROOT/'results/RUN_SUMMARY.md').open('a') as summary:
                summary.write('\n'+now()+'：启动Slingshot精确搜索v3；保留全部细胞、全部曲线顶点和原生参数；独立软件对照通过，数值/视觉验收仍待完成。\n')
            code=child.wait()
        save(LOG/'Slingshot_exit.json',dict(returncode=code,method='Slingshot',version='Slingshot_exact_search_v3'))
        successful=code==0 and (LOG/'CRC_Slingshot_completed.json').exists()
        record.update(returncode=code,finished_at=now(),elapsed_seconds=time.monotonic()-started,
                      status='fit_saved_pending_acceptance' if successful else 'fit_failed')
        save(record_file,record);save(LOG/'progress.json',record)
        with (ROOT/'results/RUN_SUMMARY.md').open('a') as summary:
            summary.write('\n'+now()+f'：Slingshot v3退出码{code}，耗时{record["elapsed_seconds"]:.1f}秒；保存状态{record["status"]}。\n')
        return 0 if successful else 1


if __name__=='__main__':
    raise SystemExit(main())
