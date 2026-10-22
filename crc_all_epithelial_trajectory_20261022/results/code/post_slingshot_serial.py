"""Wait for saved Slingshot, then run required checks and analyses one at a time."""
import json,os,subprocess,sys,time
from pipeline_common import ROOT,LOG,ENV,save_json,digest


def main():
    report=dict(status='waiting_for_formal_Slingshot',steps=[],core_analysis_accepted=False)
    record=LOG/'Slingshot_serial_postprocessing.json'
    if record.exists():raise RuntimeError('Postprocessing supervisor already recorded; inspect before restarting')
    save_json(record,report)
    while True:
        fit=json.loads((LOG/'Slingshot_local_process.json').read_text())
        if fit['status']!='running':break
        time.sleep(20)
    if fit['returncode']!=0 or fit['status']!='fit_saved_pending_acceptance':
        report['status']='fit_failed_no_postprocessing';save_json(record,report)
        raise RuntimeError('Slingshot failed; retain all logs and partial objects')
    env=os.environ.copy();env.update(TMPDIR=str(ROOT/'.runtime/tmp'),XDG_CACHE_HOME=str(ROOT/'.runtime/cache'),
        MPLCONFIGDIR=str(ROOT/'.runtime/matplotlib'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
    py=sys.executable
    steps=[('independent_Slingshot_objects',[str(ENV/'bin/Rscript'),'results/code/validate_slingshot_saved.R',str(ROOT)]),
           ('Slingshot_full_cell_figures',[py,'results/code/plot_slingshot_all.py']),
           ('Monocle2_path_dynamics',[py,'results/code/dynamics_crc.py','--method','Monocle2']),
           ('Slingshot_path_dynamics',[py,'results/code/dynamics_crc.py','--method','Slingshot']),
           ('independent_all_path_dynamics',[py,'results/code/validate_postprocess_all.py','--stage','dynamics']),
           ('independent_actual_dynamic_gene_models',[py,'results/code/validate_dynamic_model_full.py']),
           ('dynamic_root_and_all_cohort_sensitivity',[py,'results/code/dynamic_sensitivity.py']),
           ('dynamic_root_cohort_sensitivity_figures',[py,'results/code/plot_dynamic_sensitivity.py']),
           ('module_and_representative_all_path_sensitivity',[py,'results/code/module_path_dynamics.py']),
           ('MST_PCA_Harmony_and_all_cohort_sensitivity',[str(ENV/'bin/Rscript'),'results/code/topology_sensitivity.R',str(ROOT)])]
    report['status']='running_serial_postprocessing';save_json(record,report)
    for name,command in steps:
        print('START',name,flush=True);start=time.time()
        source=ROOT/command[1]
        with (LOG/(name+'_console.log')).open('x') as log:
            p=subprocess.run(['/usr/bin/time','-v','-o',str(LOG/(name+'_resources.txt'))]+command,
                cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        report['steps'].append(dict(name=name,command=command,source_sha256=digest(source),returncode=p.returncode,elapsed_seconds=time.time()-start))
        save_json(record,report);print('END',name,'exit',p.returncode,flush=True)
        if p.returncode:
            report['status']='postprocessing_failed_preserve_outputs';save_json(record,report);raise SystemExit(p.returncode)
    report['status']='saved_pending_remaining_independent_module_source_checks_and_actual_visual_review'
    save_json(record,report)


if __name__=='__main__':main()
