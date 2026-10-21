"""Assemble final acceptance only after independent checks and actual image review.

This does not fit models, replace independent validation, or approve unseen images.
"""
import datetime
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOG = ROOT / 'results/04_日志和检查'


def read(name):
    return json.loads((LOG / name).read_text())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    reports = {}
    for name in ['Monocle2_independent_validation.json',
                 'Slingshot_independent_validation.json',
                 'function_and_composition_independent_validation.json',
                 'postprocessing_independent_validation.json',
                 'CNV_export_precision_review.json',
                 'frozen_unit_hierarchy_validation.json',
                 'root_fallback_synthetic_validation.json',
                 'graph_compatibility.json',
                 'exact_projection_synthetic_validation.json',
                 'exact_entry_e2e_validation.json',
                 'final_provenance.json']:
        report = read(name)
        assert report['passed'], name
        reports[name] = digest(LOG / name)
    provenance = read('final_provenance.json')
    assert provenance['original_Monocle2_object_bytes_unchanged']
    assert not provenance['tracked_worktree_status']
    for path, expected in provenance['source_sha256'].items():
        assert digest(ROOT / path) == expected, path
    archive = read('input_archive_verification.json')
    assert archive['archive_hash_passed'] and archive['zip_crc_passed']
    source = read('frozen_source_to_work_comparison.json')
    assert source['cells'] == 150917
    for field in ['CNV_bands_exactly_equal', 'continuous_CNV_binary_exactly_equal',
                  'source_input_unchanged', 'source_and_work_count_UMI_exactly_equal']:
        assert source[field], field
    binding = read('exact_projection_software_binding.json')
    for path, expected in binding['tested_source_sha256'].items():
        assert digest(ROOT / path) == expected, path
    for name in ['input_archive_verification.json',
                 'frozen_source_to_work_comparison.json',
                 'exact_projection_software_binding.json']:
        reports[name] = digest(LOG / name)
    fit_binding = read('Slingshot_binding.json')
    assert fit_binding['original_cells_and_curve_vertices_retained']
    assert not fit_binding['fitting_parameters_changed']
    for path, expected in fit_binding['source_sha256'].items():
        assert digest(ROOT / path) == expected, path
    projection = read('Slingshot_exact_projection_audit.json')
    assert projection and all(call['input_rows_retained'] == 150917 for call in projection)
    assert all(call['curve_vertices_retained'] >= 1000 for call in projection)
    assert all(0 < call['evaluated_segment_pairs'] <= call['full_enumeration_pairs']
               for call in projection)
    for name in ['Slingshot_binding.json', 'Slingshot_exact_projection_audit.json']:
        reports[name] = digest(LOG / name)
    monocle = read('Monocle2_independent_validation.json')
    slingshot = read('Slingshot_independent_validation.json')
    for method in [monocle, slingshot]:
        for root in ['primary', 'alternative']:
            assert method[root]['cells'] == 150917
    assert not monocle['same_actual_root']
    dynamics = read('postprocessing_independent_validation.json')['dynamics']
    expected_paths = sum(method[root]['paths'] for method in [monocle, slingshot]
                         for root in ['primary', 'alternative'])
    expected_pairs = sum(monocle[root]['paths'] * slingshot[root]['paths']
                         for root in ['primary', 'alternative'])
    assert len(dynamics['paths']) == expected_paths
    assert dynamics['all_method_path_pairs_reported']
    assert dynamics['pairs'] == expected_pairs
    assert dynamics['all_library_thresholds_and_label_support_rules_verified']
    assert len(read('CNV_export_precision_review.json')['files']) == 5
    scores = read('function_and_composition_independent_validation.json')
    assert scores['scores']['all_cells_independently_rescored'] == 150917
    for path, expected in scores['verified_output_hashes'].items():
        assert digest(ROOT / path) == expected, path
    for method in ['Monocle2', 'Slingshot']:
        assert read(method + '_exit.json')['returncode'] == 0
        assert (LOG / ('CRC_' + method + '_completed.json')).exists()
    progress = read('saved_postprocess_progress.json')
    assert progress['status'] == 'numeric_checks_done_visual_review_pending'
    assert len(progress['steps']) == 6
    assert all(step['status'] == 'passed' and step['returncode'] == 0
               for step in progress['steps'])
    warnings = read('warning_review.json')
    assert warnings['status'] == 'reviewed_with_limitations'
    assert warnings['Slingshot']['status'] == 'reviewed_with_limitations'
    visual = read('visual_audit.json')
    assert visual['status'] == 'passed_all_actual_visual_reviews_bound_to_files'
    actual = {str(p.relative_to(ROOT / 'results'))
              for p in (ROOT / 'results/01_图形').rglob('*')
              if p.suffix in ['.png', '.pdf']}
    assert actual == {record['file'] for record in visual['records']}
    for record in visual['records']:
        assert record['actual_visual_inspection'] and record['status'] == 'passed'
        assert digest(ROOT / 'results' / record['file']) == record['sha256']
    for name in ['warning_review.json', 'visual_audit.json',
                 'saved_postprocess_progress.json']:
        reports[name] = digest(LOG / name)
    report = dict(
        updated_at=datetime.datetime.now().astimezone().isoformat(),
        whole_workflow_accepted=True,
        acceptance_scope='Local computational and visual delivery, with documented native warnings; Mac reception review remains separate',
        cells=150917, Monocle2=monocle, Slingshot=slingshot,
        total_method_root_paths=expected_paths, cross_method_path_pairs=expected_pairs,
        full_pool_score_max_error=scores['scores']['score_max_error'],
        png_count=visual['png_count'], pdf_count=visual['pdf_count'],
        pdf_pages_actually_reviewed=visual['pdf_page_count'],
        source_check_report_sha256=reports,
        limitations=['Native dispersion/kmeans warnings and nonmonotonic DDRTree objective reviewed; no global optimum claim',
                     'Slingshot exact segment index validated against native software fixtures; full CRC native baseline did not finish, so full-fit native equality is unavailable',
                     'Reported patient labels are not verified independent patient identities; descriptive results only',
                     'MST sensitivity uses inherited representations without refitting Harmony or curves',
                     'No new significance tests or causal/clinical conclusions'],
        pending=[])
    (LOG / 'numeric_acceptance.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
