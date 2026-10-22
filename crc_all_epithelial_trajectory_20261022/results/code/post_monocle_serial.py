"""Wait for the current fit, then execute saved-object postprocessing serially."""
import json
import os
import subprocess
import sys
import time

from pipeline_common import ROOT,LOG,ENV,save_json


def main():
    process=LOG/'Monocle2_local_process.json'
    while True:
        status=json.loads(process.read_text())
        if status['status']!='running':break
        time.sleep(2)
    if status['returncode']!=0 or status['status']!='fit_saved_pending_acceptance':
        raise RuntimeError('Formal Monocle2 did not save successfully; preserve logs and diagnose')
    env=os.environ.copy()
    env.update(TMPDIR=str(ROOT/'.runtime/tmp'),XDG_CACHE_HOME=str(ROOT/'.runtime/cache'),
               MPLCONFIGDIR=str(ROOT/'.runtime/matplotlib'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
    steps=[('independent_Monocle2_objects',[str(ENV/'bin/Rscript'),'results/code/validate_monocle_saved.R',str(ROOT)]),
           ('all_cell_function_scores',[sys.executable,'results/code/score_functions_crc.py']),
           ('main_Monocle2_figures',[sys.executable,'results/code/plot_main_all.py']),
           ('independent_scores_and_compositions',[sys.executable,'results/code/validate_postprocess_all.py','--stage','scores-and-compositions']),
           ('full_gene_bin_expression',[sys.executable,'results/code/aggregate_dynamic_expression.py']),
           ('dynamic_gene_statistics',[sys.executable,'results/code/analyze_dynamic_genes.py'])]
    report=dict(status='running_serial_postprocessing',steps=[],core_analysis_accepted=False)
    for name,command in steps:
        print('START',name,flush=True)
        started=time.time()
        with (LOG/(name+'_console.log')).open('x') as f:
            p=subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        record=dict(name=name,command=command,returncode=p.returncode,elapsed_seconds=time.time()-started)
        report['steps'].append(record)
        save_json(LOG/'Monocle2_saved_postprocessing.json',report)
        print('END',name,'exit',p.returncode,flush=True)
        if p.returncode:
            report['status']='postprocessing_failed_preserve_outputs'
            save_json(LOG/'Monocle2_saved_postprocessing.json',report)
            raise SystemExit(p.returncode)
    report['status']='saved_postprocessing_done_pending_independent_gene_checks_modules_and_visual_review'
    save_json(LOG/'Monocle2_saved_postprocessing.json',report)


if __name__=='__main__':main()
