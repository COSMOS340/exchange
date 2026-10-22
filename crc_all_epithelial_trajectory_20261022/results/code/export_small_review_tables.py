"""Keep small review tables and acceptance receipts visible on the result branch."""
import json,shutil
from pipeline_common import ROOT,TAB,LOG,digest,save_json


def main():
    numeric=json.loads((LOG/'numeric_acceptance.json').read_text())
    visual=json.loads((LOG/'visual_audit.json').read_text())
    assert numeric['passed'] and visual['passed']
    destination=ROOT/'PC_REVIEW_TABLES';destination.mkdir(exist_ok=False)
    selected=['source_counts_by_dataset.csv','source_counts_by_capture_id.csv','source_counts_by_specimen_id.csv',
        'Monocle2_roots.json','Monocle2_endpoint_root_coverage.csv','Slingshot_actual_Monocle2_reference_groups.json',
        'Slingshot_roots_and_paths.json','function_gene_coverage.csv','Monocle2_path_coverage.csv','Slingshot_path_coverage.csv',
        'topology_sensitivity.csv','CRC_only_comparison/all_four_root_comparisons.csv',
        'CRC_only_comparison/all_path_pair_overlaps.csv','source_comparisons/shared_identity_source_coverage.csv',
        'dynamic_genes/Monocle2_primary/dynamic_module_k_selection_metrics.csv',
        'dynamic_genes/Monocle2_primary/dynamic_module_bootstrap_stability.csv',
        'dynamic_genes/Monocle2_primary/dynamic_module_representative_genes.csv',
        'dynamic_genes/Monocle2_primary/dynamic_module_GO_BP_display.csv',
        'dynamic_genes/Monocle2_primary/dynamic_GO_background_coverage.csv',
        'dynamic_genes/Monocle2_primary/dynamic_model_summary.json',
        'dynamic_genes/Monocle2_alternative/dynamic_model_summary.json',
        'dynamic_genes/Monocle2_primary_all_specimen_coverage.csv',
        'dynamic_genes/Monocle2_alternative_all_specimen_coverage.csv',
        'dynamic_genes/sensitivity/all_cohort_leaveout_summary.csv',
        'dynamic_genes/method_path_sensitivity/all_method_root_path_roster.csv']
    records=[]
    for relative in selected:
        source=TAB/relative
        assert source.is_file() and source.stat().st_size<1024**2,relative
        target=destination/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        records.append(dict(path=target.relative_to(ROOT).as_posix(),canonical_result_path='results/02_表格/AllEpithelial/'+relative,
            bytes=target.stat().st_size,sha256=digest(target)))
    save_json(destination/'REVIEW_TABLE_MANIFEST.json',records)
    save_json(ROOT/'PC_ACCEPTANCE_RECEIPT.json',dict(PC_numeric_passed=True,PC_actual_visual_passed=True,
        all_cells=223736,figure_pairs=visual['formal_PNGs'],PDF_pages=visual['PDF_pages'],model_objects=numeric['model_objects'],
        canonical_numeric_acceptance_sha256=digest(LOG/'numeric_acceptance.json'),canonical_visual_audit_sha256=digest(LOG/'visual_audit.json'),
        Monocle2_actual_fit_code_commit='bdcc8df4fdf927b9d2db11ad349c1e60da76b8d8',
        Slingshot_actual_fit_code_commit='adeb6886ff58a7779dc36ac8e961eb6937112691',
        final_delivery_code_commit='Recorded after clean commit by RETURN_MANIFEST.json and exact Release tag target',
        full_output_files_and_hashes='RESULT_MANIFEST.json in new result Release',
        local_transport_and_remote_publication_verification='Separate receipts after packing/uploading',
        Mac_receiving_acceptance='pending; do not automatically replace the formal step14 entry'))
    (destination/'README.md').write_text('# PC小型核验表\n\n这些文件逐字复制自新结果包的规范源表，路径和SHA见REVIEW_TABLE_MANIFEST.json。完整逐细胞/全基因/分箱表、图形、RDS对象、环境及日志在新的结果Release。PC数值/视觉检查通过不代表Mac已接收。\n')
    print('Small result-branch review tables and PC acceptance receipt exported.',flush=True)


if __name__=='__main__':main()
