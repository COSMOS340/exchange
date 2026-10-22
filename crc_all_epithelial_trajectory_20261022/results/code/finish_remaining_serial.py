"""Run remaining independent checks after the live serial supervisor has ended."""
import json,os,subprocess,sys,time
from pipeline_common import ROOT,LOG,ENV,save_json,digest


def main():
    destination=LOG/'manual_remaining_steps.json'
    if destination.exists():raise RuntimeError('A remaining-steps supervisor is already recorded')
    report=dict(status='waiting_for_previous_serial_supervisor',steps=[],scientific_acceptance=False)
    save_json(destination,report)
    waiting={'waiting_for_formal_Slingshot','running_serial_postprocessing'}
    while True:
        previous=json.loads((LOG/'Slingshot_serial_postprocessing.json').read_text())
        if previous['status'] not in waiting:break
        time.sleep(20)
    expected='saved_pending_remaining_independent_module_source_checks_and_actual_visual_review'
    if previous['status']!=expected:
        report['status']='previous_stage_failed_do_not_continue';save_json(destination,report)
        raise RuntimeError(previous['status'])
    env=os.environ.copy();env.update(TMPDIR=str(ROOT/'.runtime/tmp'),XDG_CACHE_HOME=str(ROOT/'.runtime/cache'),
        MPLCONFIGDIR=str(ROOT/'.runtime/matplotlib'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
    steps=[('independent_all_cell_module_values_and_paths',[sys.executable,'results/code/validate_module_activity_paths.py']),
        ('independent_all_cohort_dynamic_sensitivity',[sys.executable,'results/code/validate_dynamic_sensitivity.py']),
        ('remaining_root_cohort_module_figures',[sys.executable,'results/code/plot_dynamic_sensitivity.py']),
        ('native_warning_evidence',[sys.executable,'results/code/collect_warning_evidence.py']),
        ('all_remaining_PDF_rendering',['/home/anran/hehe/.crc-visual-audit-env/bin/python','results/code/render_pdf_audit.py'])]
    report['status']='running_remaining_serial_checks';save_json(destination,report)
    for name,command in steps:
        start=time.time();print('START',name,flush=True)
        with (LOG/(name+'_console.log')).open('x') as log:
            process=subprocess.run(['/usr/bin/time','-v','-o',str(LOG/(name+'_resources.txt'))]+command,
                cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        report['steps'].append(dict(name=name,command=command,source_sha256=digest(ROOT/command[1]),returncode=process.returncode,elapsed_seconds=time.time()-start))
        save_json(destination,report);print('END',name,'exit',process.returncode,flush=True)
        if process.returncode:
            report['status']='failed_preserve_outputs';save_json(destination,report);raise SystemExit(process.returncode)
    report['status']='independent_remaining_checks_done_PENDING_ACTUAL_VISUAL_REVIEW_AND_WARNING_REVIEW'
    save_json(destination,report)


if __name__=='__main__':main()
