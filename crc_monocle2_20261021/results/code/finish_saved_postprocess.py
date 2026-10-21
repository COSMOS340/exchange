"""Wait for this checkout's running fit, then serially check its saved outputs.

Never starts or restarts a model fit, publishes results, or approves visual QA.
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

ROOT = Path(__file__).resolve().parents[2]
LOG = ROOT / 'results/04_日志和检查'


def stamp():
    return datetime.datetime.now().astimezone().isoformat()


def save(record):
    temporary = LOG / 'saved_postprocess_progress.tmp'
    temporary.write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n')
    temporary.replace(LOG / 'saved_postprocess_progress.json')


def main():
    with (LOG/'saved_postprocess.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (LOG/'saved_postprocess_progress.json').exists():
            raise RuntimeError('Prior postprocess record exists; inspect and version any retry')
        record = dict(started_at=stamp(), phase='waiting_for_Slingshot',
                      status='running', steps=[], whole_analysis_accepted=False)
        save(record)
        print(stamp(), record['phase'], flush=True)
        while True:
            fit = json.loads((LOG/'Slingshot_local_process.json').read_text())
            if fit['status'] != 'running':
                break
            time.sleep(15)
        if fit['returncode'] != 0 or not (LOG/'CRC_Slingshot_completed.json').exists():
            record.update(status='blocked_by_fit_failure', finished_at=stamp(), fit=fit)
            save(record)
            raise RuntimeError('Slingshot did not save a complete model; preserve diagnostics')
        # The model supervisor holds fit.lock until its child has exited.
        with (LOG/'fit.lock').open('a') as fit_lock:
            fcntl.flock(fit_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            code = ROOT/'results/code'
            rscript = ROOT/'.venv/bin/Rscript'
            commands = [
                ('validate_slingshot', [str(rscript), str(code/'validate_slingshot_saved.R'), str(ROOT)]),
                ('MST_sensitivity', [str(rscript), str(code/'topology_sensitivity.R'), str(ROOT)]),
                ('plot_slingshot', [sys.executable, '-u', str(code/'plot_slingshot_crc.py')]),
                ('Slingshot_dynamics', [sys.executable, '-u', str(code/'dynamics_crc.py'), '--method', 'Slingshot']),
                ('validate_all_dynamics', [sys.executable, '-u', str(code/'validate_postprocess_crc.py'), '--stage', 'dynamics']),
                ('render_all_PDFs', ['/home/anran/hehe/.crc-visual-audit-env/bin/python', str(code/'render_pdf_audit.py')]),
            ]
            env = os.environ.copy()
            for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
                env[key] = '1'
            env['MPLCONFIGDIR'] = str(ROOT.parent.parent/'.crc-cache/matplotlib')
            for name, command in commands:
                script = next(Path(arg) for arg in command if arg.endswith(('.R','.py')))
                step = dict(name=name, command=command, cwd=str(ROOT), started_at=stamp(),
                            code_sha256=hashlib.sha256(script.read_bytes()).hexdigest(), status='running')
                record['steps'].append(step); record['phase'] = name; save(record)
                print(stamp(), name, flush=True)
                with (LOG/f'{name}_serial_console.log').open('x') as output:
                    result = subprocess.run(command, cwd=ROOT, env=env, stdout=output, stderr=subprocess.STDOUT)
                step.update(returncode=result.returncode, finished_at=stamp(),
                            status='passed' if result.returncode == 0 else 'failed')
                save(record)
                if result.returncode:
                    record.update(status='postprocessing_failed', finished_at=stamp()); save(record)
                    raise RuntimeError(f'{name} failed; inspect its retained console before retry')
        record.update(status='numeric_checks_done_visual_review_pending',
                      phase='actual_PNG_and_every_PDF_page_review', finished_at=stamp())
        save(record)
        print(stamp(), record['status'], flush=True)


if __name__ == '__main__':
    main()
