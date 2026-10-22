"""Require actual per-file visual observations and all numerical acceptance evidence."""
import json
from pathlib import Path
import pandas as pd
from pipeline_common import ROOT,TAB,OBJ,LOG,digest,save_json


def figure_records():
    records={}
    for path in sorted(LOG.glob('*manifest.json')):
        value=json.loads(path.read_text())
        if not isinstance(value,list):continue
        for row in value:
            if not isinstance(row,dict) or not {'name','png','pdf'}<=set(row):continue
            name=row['name']
            if name in records:
                assert records[name]['png']==row['png'] and records[name]['pdf']==row['pdf']
            records.setdefault(name,{}).update(row)
    return records


def visual_acceptance():
    records=figure_records()
    observations=json.loads((LOG/'actual_visual_observations.json').read_text())
    pages=json.loads((LOG/'pdf_page_manifest.json').read_text())
    pngs={p.relative_to(ROOT/'results').as_posix() for p in (ROOT/'results/01_图形').rglob('*.png')}
    pdfs={p.relative_to(ROOT/'results').as_posix() for p in (ROOT/'results/01_图形').rglob('*.pdf')}
    assert pngs=={r['png'] for r in records.values()}
    assert pdfs=={r['pdf'] for r in records.values()}=={p['pdf'] for p in pages}
    verified=[]
    for name,row in sorted(records.items()):
        observation=observations[name]
        assert observation['png_actually_viewed'] is True and observation['pdf_actually_viewed'] is True,name
        assert len(observation['observation'])>=25,name
        for kind in ['png','pdf']:
            assert digest(ROOT/'results'/row[kind])==observation[kind+'_sha256'],name+' '+kind
        own_pages=[p for p in pages if p['pdf']==row['pdf']]
        assert own_pages and len(own_pages)==own_pages[0]['total_pages']
        assert set(observation['pdf_pages_viewed'])=={p['page'] for p in own_pages}
        for page in own_pages:
            assert page['source_pdf_sha256']==observation['pdf_sha256']
            assert page['rendered_png_sha256']==digest(ROOT/'results'/page['rendered_png'])
        for source in row.get('sources',[]):assert (ROOT/'results'/source).is_file(),source
        verified.append(dict(name=name,png=row['png'],pdf=row['pdf'],pages=len(own_pages),
            png_sha256=observation['png_sha256'],pdf_sha256=observation['pdf_sha256'],observation=observation['observation']))
    assert any(r['name'].startswith('12_Slingshot') for r in verified),'Slingshot figure review still missing'
    assert any(r['name'].startswith('50_') for r in verified),'Module method sensitivity figure review still missing'
    assert any(r['name'].startswith('34_') for r in verified),'All-cohort module sensitivity figure review still missing'
    assert any(r['name'].endswith('supported_labels') for r in verified),'Path activity figure review still missing'
    # Synchronize status only after explicit observations for ALL current files
    # have passed the checks above. Rendering/file presence is not an attestation.
    for path in sorted(LOG.glob('*manifest.json')):
        value=json.loads(path.read_text())
        if not isinstance(value,list):continue
        changed=False
        for row in value:
            if not isinstance(row,dict) or not {'name','png','pdf'}<=set(row):continue
            observation=observations[row['name']]
            row.update(visual_review='passed_actual_opening_and_SHA_verified',
                actual_visual_observation=observation['observation'],
                png_sha256=observation['png_sha256'],pdf_sha256=observation['pdf_sha256'])
            if 'Pending' in str(row.get('finding','')) or not row.get('finding'):row['finding']=observation['observation']
            changed=True
        if changed:save_json(path,value)
    for page in pages:page['visual_review']='passed_actual_page_opening_and_source_SHA_verified'
    save_json(LOG/'pdf_page_manifest.json',pages)
    save_json(LOG/'visual_audit.json',dict(passed=True,formal_PNGs=len(pngs),formal_PDFs=len(pdfs),PDF_pages=len(pages),
        actual_per_file_observations_not_inferred_from_rendering=True,figures=verified))
    return records,observations


def main():
    checks=['input_archive_validation.json','preparation_validation.json','python_unit_test_result.json',
        'graph_compatibility_validation.json','root_fallback_synthetic_validation.json','endpoint_root_synthetic_validation.json',
        'exact_projection_synthetic_validation.json','exact_entry_e2e_validation.json',
        'dynamic_model_software_tests.json','module_clustering_software_tests.json',
        'Monocle2_independent_validation.json','Slingshot_independent_validation.json',
        'function_and_composition_independent_validation.json','postprocessing_independent_validation.json',
        'dynamic_expression_independent_validation.json','dynamic_model_full_independent_validation.json',
        'source_profiles_and_module_independent_validation.json','GO_independent_validation.json',
        'CRC_only_comparison_validation.json','CRC_comparison_independent_validation.json','dynamic_sensitivity_independent_validation.json',
        'module_activity_path_independent_validation.json','all_specimen_dynamic_coverage.json']
    evidence=[]
    for name in checks:
        p=LOG/name;value=json.loads(p.read_text())
        if name=='input_archive_validation.json':
            assert value['archive_verified'] is True and value['archive_CRC_verified'] is True,name
        else:assert value.get('passed') is True,name
        evidence.append(dict(path='04_日志和检查/'+name,sha256=digest(p)))
    for method in ['Monocle2','Slingshot']:
        process=json.loads((LOG/(method+'_local_process.json')).read_text())
        assert process['returncode']==0 and process['status']=='fit_saved_pending_acceptance'
    warning=json.loads((LOG/'warning_and_convergence_review.json').read_text())
    assert warning['reviewed'] is True and warning['limitations_written_in_methods'] is True
    references=json.loads((TAB/'Slingshot_actual_Monocle2_reference_groups.json').read_text())
    roots=json.loads((TAB/'Monocle2_roots.json').read_text())
    for root in ['primary','alternative']:
        assert references[root]==roots[root]['actual_root_group']
        cells=pd.read_csv(TAB/f'Monocle2_{root}_cells.csv.gz')
        assert len(cells)==223736 and cells.cell_id.is_unique
    expected=['Monocle2_checkpoint.rds','Monocle2_primary.rds','Monocle2_alternative.rds',
        'Slingshot_primary.rds','Slingshot_alternative.rds','Slingshot_primary_UMAP.rds','Slingshot_alternative_UMAP.rds']
    # Use the actual checkpoint name recorded by the native engine if different.
    if not (OBJ/expected[0]).exists():
        candidates=list(OBJ.glob('*checkpoint*.rds'));assert len(candidates)==1;expected[0]=candidates[0].name
    objects=[]
    for name in expected:
        p=OBJ/name;assert p.is_file() and p.stat().st_size>1024**2
        objects.append(dict(path='03_对象/AllEpithelial/'+name,bytes=p.stat().st_size,sha256=digest(p)))
    save_json(LOG/'model_object_SHA256.json',objects)
    records,observations=visual_acceptance()
    save_json(LOG/'numeric_acceptance.json',dict(passed=True,cells=223736,
        two_class_counts={'Normal_colorectal':72819,'Primary_CRC':150917},
        three_class_counts={'Normal':64654,'Adjacent':8165,'CRC':150917},
        captures=242,specimens=229,reported_subject_labels=149,verified_independent_patients=False,
        checks=evidence,model_objects=objects,formal_figure_pairs=len(records),
        interpretation='Numerical consistency and scoped PC deliverables accepted with explicit native warnings; no proof of clinical time, progression, causality or independent biological validation',
        Mac_receiving_acceptance='pending; do not update the formal step14 entry automatically'))
    print('All required numerical, object-hash and actual visual evidence passed.',flush=True)


if __name__=='__main__':main()
