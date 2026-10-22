"""Inventory native warnings and diagnostics; never auto-approve a fit."""
import json,re
from collections import Counter
from pipeline_common import ROOT,LOG,save_json,digest


def main():
    report=dict(reviewed=False,limitations_written_in_methods=False,methods={})
    for method in ['Monocle2','Slingshot']:
        process=json.loads((LOG/f'{method}_local_process.json').read_text())
        assert process['status']=='fit_saved_pending_acceptance' and process['returncode']==0
        path=LOG/f'{method}_console.log';lines=path.read_text(errors='replace').splitlines()
        warnings=[]
        for i,line in enumerate(lines):
            if re.search(r'Warning|did not converge|steps exceeded|NaNs produced',line):
                warnings.append(dict(line=i+1,text=line,following=lines[i+1:i+3]))
        resource=(LOG/f'{method}_resource_usage.txt').read_text()
        report['methods'][method]=dict(console_sha256=digest(path),native_warnings=warnings,
            GNU_time_resource_text=resource,process=process,
            independent_saved_object_validation_sha256=digest(LOG/f'{method}_independent_validation.json'))
    report['Monocle2_interpretation']=dict(
        dispersion_GLMS_and_parametric_fit_warning='Retained. Finite saved dispersions do not prove every GLM/parametric fit converged.',
        dense_conversion='Native 1998-ordering-gene expression panel, about3.3GiB; not a cell N-by-N matrix.',
        Quick_TRANSFER_limit='Retained native kmeans warning; no guarantee its transfer stage fully converged.',
        DDRTree_objective=json.loads((LOG/'Monocle2_independent_validation.json').read_text())['objective'],
        global_optimum_or_all_native_optimizers_converged_claimed=False)
    report['Slingshot_interpretation']=dict(max_iterations=30,threshold=.001,
        native_termination='abs(old_distance-new_distance) <= threshold*old_distance OR maxit reached',
        native_source_sha256=digest(LOG/'native_Slingshot_getCurves_source.txt'),
        iteration_history_or_hasConverged_flag_persisted=False,
        limitation='Native function returns curves/PT/weights without iteration count or hasConverged in the object; success alone cannot identify which stopping criterion ended the loop.',
        approximation_points=150,all_original_cells_and_initial_projection_segments_retained=True,
        unaccelerated_full_biological_baseline_repeated=False)
    report['all_native_warning_text_preserved']=True
    report['required_manual_action']='Read raw warning contexts and independent object/parameter/coverage diagnostics; describe their limits in methods, then explicitly record the review decision. This collector does not set reviewed=true.'
    save_json(LOG/'warning_evidence_pending_review.json',report)
    print('Native warning/convergence evidence collected; manual review still required.',flush=True)


if __name__=='__main__':main()
