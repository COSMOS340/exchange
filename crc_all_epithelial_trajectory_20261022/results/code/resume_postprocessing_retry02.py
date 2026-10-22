"""Resume only the failed probability validator and previously unstarted steps."""
import os,subprocess,sys,time
from pipeline_common import ROOT,LOG,ENV,save_json,digest


def main():
    report=dict(status='running_retry02_after_explicit_R_boolean_import_fix',steps=[],scientific_acceptance=False)
    destination=LOG/'Slingshot_postprocessing_retry02.json'
    if destination.exists():raise RuntimeError('Retry02 already recorded; inspect before repeating')
    env=os.environ.copy();env.update(TMPDIR=str(ROOT/'.runtime/tmp'),XDG_CACHE_HOME=str(ROOT/'.runtime/cache'),
        MPLCONFIGDIR=str(ROOT/'.runtime/matplotlib'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
    py=sys.executable
    steps=[('independent_actual_dynamic_gene_models_retry02',[py,'results/code/validate_dynamic_model_full.py']),
        ('dynamic_root_and_all_cohort_sensitivity',[py,'results/code/dynamic_sensitivity.py']),
        ('module_and_representative_all_path_sensitivity',[py,'results/code/module_path_dynamics.py']),
        ('MST_PCA_Harmony_and_all_cohort_sensitivity',[str(ENV/'bin/Rscript'),'results/code/topology_sensitivity.R',str(ROOT)]),
        ('independent_all_cell_module_values_and_paths',[py,'results/code/validate_module_activity_paths.py']),
        ('independent_all_cohort_dynamic_sensitivity',[py,'results/code/validate_dynamic_sensitivity.py']),
        ('remaining_root_cohort_module_figures',[py,'results/code/plot_dynamic_sensitivity.py']),
        ('native_warning_evidence',[py,'results/code/collect_warning_evidence.py']),
        ('all_remaining_PDF_rendering',['/home/anran/hehe/.crc-visual-audit-env/bin/python','results/code/render_pdf_audit.py'])]
    save_json(destination,report)
    for name,command in steps:
        print('START',name,flush=True);start=time.time();source_sha=digest(ROOT/command[1])
        with (LOG/(name+'_console.log')).open('x') as log:
            process=subprocess.run(['/usr/bin/time','-v','-o',str(LOG/(name+'_resources.txt'))]+command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        report['steps'].append(dict(name=name,command=command,source_sha256=source_sha,returncode=process.returncode,elapsed_seconds=time.time()-start))
        save_json(destination,report);print('END',name,'exit',process.returncode,flush=True)
        if process.returncode:
            report['status']='retry02_failed_preserve_outputs';save_json(destination,report);raise SystemExit(process.returncode)
    report['status']='all_remaining_computations_done_PENDING_ACTUAL_VISUAL_AND_WARNING_REVIEW';save_json(destination,report)


if __name__=='__main__':main()
