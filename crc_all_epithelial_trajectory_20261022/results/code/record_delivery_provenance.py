"""Bind delivered code, frozen task/plan and actual per-stage execution receipts."""
import datetime,json,subprocess
from pipeline_common import ROOT,LOG,TAB,save_json,digest


def main():
    for method in ['Monocle2','Slingshot']:
        process=json.loads((LOG/f'{method}_local_process.json').read_text())
        assert process['returncode']==0 and process['status']=='fit_saved_pending_acceptance'
    code=[]
    for p in sorted((ROOT/'results/code').rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ['.pyc','.pyo']:
            code.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=digest(p)))
    configuration=[]
    for name in ['INPUT_MANIFEST.json','WORK_MANIFEST.json','ANALYSIS_PLAN.json','config.json','TASK_SPEC.md',
                 'SOURCE_TISSUE_AUDIT.md','prior_specimen_source_audit.csv','REFERENCE_FIGURE_REQUIREMENTS.md','AllEpithelium_PC.py']:
        p=ROOT/name
        if p.exists():configuration.append(dict(path=name,bytes=p.stat().st_size,sha256=digest(p)))
    configuration.append(dict(path='work/fit_generated.R',bytes=(ROOT/'work/fit_generated.R').stat().st_size,sha256=digest(ROOT/'work/fit_generated.R')))
    records=[]
    names=['source_binding.json','engine_derivation.json','runtime_environment.json','Monocle2_binding.json','Slingshot_binding.json',
        'Monocle2_local_process.json','Slingshot_local_process.json','Monocle2_saved_postprocessing.json',
        'Monocle2_extension_process.json','Slingshot_serial_postprocessing.json','Slingshot_postprocessing_retry02.json','warning_and_convergence_review.json',
        'model_object_SHA256.json','numeric_acceptance.json','visual_audit.json','manual_remaining_steps.json']
    for name in names:
        p=LOG/name
        if p.exists():records.append(dict(path='results/04_日志和检查/'+name,sha256=digest(p)))
    save_json(LOG/'DELIVERY_PROVENANCE.json',dict(
        created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),task_repository='https://github.com/COSMOS340/exchange',
        task_commit='5cfd63241160058869ccc297bf4222b439813487',branch='pc-all-epithelial-results-20261022',
        document_generation_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        actual_Monocle_fit_commit='bdcc8df4fdf927b9d2db11ad349c1e60da76b8d8',actual_Slingshot_fit_commit='adeb6886ff58a7779dc36ac8e961eb6937112691',
        running_post_Slingshot_supervisor_compiled_commit='97f240a',
        individual_postprocessing_commands_and_source_hashes_in_receipts=True,
        delivery_commit='Exact clean final code+docs commit is pinned by RETURN_MANIFEST.json and the new Release tag; no circular self-hash is fabricated here',
        cells=223736,genes=15635,raw_integer_nnz=497539060,ordering_genes=1998,
        original_input_release='crc-monocle2-input-v1-20261021',input_archive_bytes=1181830179,
        input_archive_sha256='ce89b556ea4c167c2f164ed6757ec54094d9c712ff81d4db9daf4221e013d724',
        source_audit_join_unmatched=0,old_CRC_models_preserved=True,
        code_files=code,frozen_configuration=configuration,execution_and_acceptance_receipts=records,
        full_file_inventory='RESULT_MANIFEST.json (all result files) and SIDECAR_MANIFEST.json (all bound root/work sidecars) in result Release',
        all_statistics_are_conditional_or_descriptive=True,Mac_receiving_acceptance='pending'))
    print('Delivery code/plan/input and actual execution provenance bound.',flush=True)


if __name__=='__main__':main()
