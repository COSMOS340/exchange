"""Independent full-pool arithmetic and NumPy unit rollups of saved postprocessing."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/02_表格/AllEpithelial"
LOG = ROOT / "results/04_日志和检查"
SETS = json.loads((ROOT / "config.json").read_text())["function_sets"]
FEATURES = ["cnv_abs_log2", "CNV_Low_fraction", "CNV_Intermediate_fraction", "CNV_High_fraction",
            *SETS, "cycling_fraction", "log1p_UMI"]


def equal(actual, expected, tolerance=1e-8):
    a, b = np.asarray(actual, float), np.asarray(expected, float)
    assert a.shape == b.shape and np.allclose(a, b, atol=tolerance, rtol=1e-10, equal_nan=True)
    valid = np.isfinite(a) & np.isfinite(b)
    return float(np.abs(a[valid]-b[valid]).max()) if valid.any() else 0.


def aligned(path, metadata):
    data = pd.read_csv(path)
    assert data.cell_id.is_unique and set(data.cell_id) == set(metadata.cell_id)
    return data.set_index("cell_id").loc[metadata.cell_id].reset_index()


def aggregate(frame, keys, columns):
    """Factorization/addition implementation, independent of authoring groupby.mean."""
    codes, levels = pd.factorize(pd.MultiIndex.from_frame(frame[keys]), sort=False)
    assert (codes >= 0).all()
    values = frame[columns].to_numpy(dtype=float)
    valid = np.isfinite(values)
    sums = np.zeros((len(levels), len(columns)))
    counts = np.zeros_like(sums, dtype=np.int64)
    np.add.at(sums, codes, np.where(valid, values, 0))
    np.add.at(counts, codes, valid.astype(np.int64))
    means = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=counts>0)
    result = levels.to_frame(index=False)
    # pandas.factorize may discard MultiIndex level names; restore input keys.
    result.columns = keys
    result[columns] = means
    return result, np.bincount(codes, minlength=len(levels)), codes


def compare_units(expected, path, keys, columns):
    saved = pd.read_csv(path)
    for key in keys:
        if pd.api.types.is_string_dtype(expected[key].dtype):
            saved[key] = saved[key].astype(str)
    assert not saved.duplicated(keys).any() and len(saved) == len(expected)
    index = pd.MultiIndex.from_frame(expected[keys]) if len(keys)>1 else pd.Index(expected[keys[0]],name=keys[0])
    saved = saved.set_index(keys).reindex(index)
    return equal(saved[columns], expected[columns])


def validate_scores(metadata, values):
    parameters = pd.read_csv(TAB / "function_gene_standardization.csv").set_index("gene")
    members = pd.read_csv(TAB / "function_gene_membership.csv")
    coverage = pd.read_csv(TAB / "function_gene_coverage.csv").set_index("pathway")
    genes = pd.Index(pd.read_csv(ROOT / "work/AllEpithelial/genes.csv").gene)
    panel = parameters.index
    indices = genes.get_indexer(panel)
    assert (indices >= 0).all()
    mean, sd = parameters.mean_log1p_CP10K.to_numpy(), parameters.population_sd.to_numpy()
    sums, squares = np.zeros(len(panel)), np.zeros(len(panel))
    masks = {name: panel.isin(members.loc[members.pathway.eq(name) & members.used_for_score, "gene"]) for name in SETS}
    maximum = 0.; tested = 0
    with h5py.File(ROOT / "work/AllEpithelial/counts.h5") as h:
        shape = tuple(map(int, h["shape"][:]))
        pointers = h["indptr"][:]
        for lo in range(0, len(metadata), 2000):
            hi = min(lo+2000, len(metadata)); a, b = int(pointers[lo]), int(pointers[hi])
            # Native genes-by-cells CSC, unlike the scoring script's CSR transpose.
            counts = sp.csc_matrix((h["data"][a:b], h["indices"][a:b], pointers[lo:hi+1]-a),
                                   shape=(shape[0], hi-lo))
            assert np.array_equal(np.asarray(counts.sum(axis=0)).ravel(), metadata.common_gene_UMI.iloc[lo:hi])
            dense = np.log1p(counts[indices, :].toarray().T * (1e4/metadata.nCount_input.iloc[lo:hi].to_numpy())[:, None])
            sums += dense.sum(axis=0); squares += (dense*dense).sum(axis=0)
            for name in SETS:
                if coverage.loc[name, "assessable"]:
                    mask = masks[name]
                    assert mask.sum() >= 10
                    reconstructed = ((dense[:,mask]-mean[mask])/sd[mask]).mean(axis=1)
                    maximum = max(maximum, equal(values[name].iloc[lo:hi], reconstructed))
                else:
                    assert values[name].iloc[lo:hi].isna().all()
            tested += hi-lo
    recomputed_mean = sums/tested
    recomputed_sd = np.sqrt(np.maximum(0, squares/tested-recomputed_mean**2))
    mean_error, sd_error = equal(mean, recomputed_mean), equal(sd, recomputed_sd)
    assert tested == 223736
    equal(values.log1p_UMI, np.log1p(metadata.nCount_input))
    equal(values.cnv_abs_log2, metadata.cnv_abs_log2, tolerance=1e-14)
    for band in ["Low", "Intermediate", "High"]:
        equal(values["CNV_"+band+"_fraction"], metadata.CNV_relative_band.eq(band).astype(int), tolerance=0)
    return dict(all_cells_independently_rescored=tested, score_max_error=maximum,
                population_mean_max_error=mean_error, population_sd_max_error=sd_error,
                denominator="nCount_input: original native total measured UMI")


def validate_compositions(metadata):
    primary = aligned(TAB / "Monocle2_primary_cells.csv.gz", metadata)
    cases = [("Monocle2_State", metadata.assign(State=primary.State.astype(str)), "State")]
    for orientation in ["primary", "alternative"]:
        paths = aligned(TAB / f"Monocle2_{orientation}_path_cells.csv.gz", metadata)
        pieces = []
        for field in [c for c in paths if c.startswith("PT_")]:
            number = field.split("_")[1]
            keep = paths[field].notna() & paths["W_"+number].ge(.5)
            piece = metadata.loc[keep].copy(); piece["path"] = "Path"+number; pieces.append(piece)
        cases.append((f"Monocle2_{orientation}_paths", pd.concat(pieces, ignore_index=True), "path"))
    records=[]
    bands=["Low","Intermediate","High"]
    for prefix, frame, key in cases:
        frame=frame.copy()
        for band in bands:frame[band]=frame.CNV_relative_band.eq(band).astype(float)
        pooled, counts, _ = aggregate(frame,[key],bands)
        captures, _, _ = aggregate(frame,[key,"dataset","specimen_id","capture_id"],bands)
        unit, _, _ = aggregate(captures,[key,"dataset","specimen_id"],bands)
        specimens, _, _ = aggregate(unit,[key],bands)
        pooled_error=compare_units(pooled,TAB/f"{prefix}_CNV_pooled_fraction.csv",[key],bands)
        specimen_error=compare_units(specimens,TAB/f"{prefix}_CNV_specimen_equal_fraction.csv",[key],bands)
        records.append(dict(table=prefix,pooled_error=pooled_error,equal_specimen_error=specimen_error))
    return records


def validate_dynamics(metadata, values, methods=('Monocle2','Slingshot')):
    folder=TAB/'dynamics'
    source={}
    for method in methods:
        for root in ['primary','alternative']:
            suffix='path_cells' if method=='Monocle2' else 'cells'
            source[method,root]=aligned(TAB/f'{method}_{root}_{suffix}.csv.gz',metadata)
    records=[]
    expected_paths=set()
    for (method,root),paths in source.items():
        for field in [x for x in paths if x.startswith('PT_')]:
            n=field.split('_')[1]
            valid=paths[field].notna() & paths['W_'+n].ge(.5)
            if valid.any() and np.ptp(paths.loc[valid,field].to_numpy())>0:
                expected_paths.add(f'{method}_{root}_Path{n}')
    for filename in sorted(folder.glob('*_cell_bins.csv.gz')):
        stem=filename.name.removesuffix('_cell_bins.csv.gz')
        method,root,path=stem.split('_');number=path.removeprefix('Path')
        if method not in methods:
            continue
        paths=source[method,root]
        pt=paths['PT_'+number].to_numpy();weight=paths['W_'+number].to_numpy()
        valid=np.isfinite(pt)&(weight>=.5)
        assignment=pd.read_csv(filename)
        assert np.array_equal(assignment.cell_id.to_numpy(),metadata.cell_id.to_numpy()[valid])
        scaled=(pt[valid]-pt[valid].min())/np.ptp(pt[valid])
        equal(assignment.raw_pseudotime,pt[valid]);equal(assignment.scaled_pseudotime,scaled)
        equal(assignment.weight,weight[valid])
        bins=np.minimum(np.floor(12*scaled).astype(int),11)
        assert np.array_equal(assignment.bin,bins)
        boundaries=pd.read_csv(folder/f'{stem}_boundaries.csv')
        equal(boundaries.scaled_left,np.arange(12)/12)
        equal(boundaries.scaled_right,np.arange(1,13)/12)
        equal(boundaries.raw_left,pt[valid].min()+np.arange(12)/12*np.ptp(pt[valid]))
        equal(boundaries.raw_right,pt[valid].min()+np.arange(1,13)/12*np.ptp(pt[valid]))
        assert boundaries.right_inclusive.tolist()==[False]*11+[True]
        keys=['dataset','subject_key','specimen_id','capture_id','reported_patient_link_supported','bin']
        frame=metadata.loc[valid,keys[:-1]].reset_index(drop=True).copy()
        frame['bin']=bins;frame[FEATURES]=values.loc[valid,FEATURES].to_numpy()
        libraries,cells,_=aggregate(frame,keys,FEATURES)
        libraries['cells']=cells;libraries['eligible']=cells>=10
        library_error=compare_units(libraries,folder/f'{stem}_library_bins.csv.gz',keys,FEATURES+['cells','eligible'])
        eligible=libraries[libraries.eligible]
        spec_keys=['dataset','subject_key','specimen_id','reported_patient_link_supported','bin']
        specimens,n_libraries,codes=aggregate(eligible,spec_keys,FEATURES)
        specimens['libraries']=n_libraries
        spec_cells=np.zeros(len(specimens));np.add.at(spec_cells,codes,eligible.cells.to_numpy())
        specimens['cells']=spec_cells
        specimen_error=compare_units(specimens,folder/f'{stem}_specimen_bins.csv.gz',spec_keys,FEATURES+['cells','libraries'])
        unsupported=specimens[~specimens.reported_patient_link_supported]
        compare_units(unsupported,folder/f'{stem}_unsupported_specimens.csv.gz',spec_keys,FEATURES+['cells','libraries'])
        supported=specimens[specimens.reported_patient_link_supported]
        label_keys=['dataset','subject_key','bin']
        labels,n_specimens,codes=aggregate(supported,label_keys,FEATURES)
        labels['specimens']=n_specimens
        for column in ['cells','libraries']:
            totals=np.zeros(len(labels));np.add.at(totals,codes,supported[column].to_numpy());labels[column]=totals
        label_error=compare_units(labels,folder/f'{stem}_reported_label_bins.csv.gz',label_keys,FEATURES+['cells','libraries','specimens'])
        bin_coverage=pd.read_csv(folder/f'{stem}_bin_coverage.csv').set_index('bin')
        assert bin_coverage.index.tolist()==list(range(12))
        for bin_id in range(12):
            raw=frame[frame.bin.eq(bin_id)];cap=eligible[eligible.bin.eq(bin_id)]
            spec=specimens[specimens.bin.eq(bin_id)];lab=labels[labels.bin.eq(bin_id)]
            proxy=unsupported[unsupported.bin.eq(bin_id)]
            expected=dict(member_cells=len(raw),member_captures=raw.capture_id.nunique(),
                member_specimens=raw.specimen_id.nunique(),member_datasets=raw.dataset.nunique(),
                eligible_cells=cap.cells.sum(),eligible_captures=len(cap),eligible_specimens=len(spec),eligible_datasets=spec.dataset.nunique(),
                supported_reported_labels=len(lab),supported_label_cells=lab.cells.sum(),supported_label_captures=lab.libraries.sum(),
                supported_label_specimens=lab.specimens.sum(),supported_label_datasets=lab.dataset.nunique(),
                unsupported_proxy_specimens=len(proxy),unsupported_proxy_cells=proxy.cells.sum())
            equal(bin_coverage.loc[bin_id,list(expected)],list(expected.values()))
        changes=pd.read_csv(TAB/f'{method}_within_reported_label_changes.csv.gz')
        changes=changes[changes.root.eq(root)&changes.path.eq(int(number))]
        checked_changes=0
        for (dataset,subject),group in labels.groupby(['dataset','subject_key']):
            for feature in FEATURES:
                measured=group[group[feature].notna()].sort_values('bin')
                if len(measured)<2:continue
                row=changes[changes.dataset.eq(dataset)&changes.subject_key.eq(subject)&changes.feature.eq(feature)]
                assert len(row)==1 and row.observed_bins.iloc[0]==len(measured)
                assert row.first_bin.iloc[0]==measured.bin.iloc[0] and row.last_bin.iloc[0]==measured.bin.iloc[-1]
                equal(row[['first_value','last_value','difference']],np.array([[measured[feature].iloc[0],measured[feature].iloc[-1],measured[feature].iloc[-1]-measured[feature].iloc[0]]]))
                checked_changes+=1
        assert checked_changes==len(changes)
        summary=pd.read_csv(folder/f'{stem}_summary.csv').set_index(['bin','feature'])
        summary_error=0.;support_count=0
        for bin_id in range(12):
            for feature in FEATURES:
                x=labels.loc[labels.bin.eq(bin_id),feature].dropna().to_numpy()
                row=summary.loc[(bin_id,feature)]
                assert int(row.reported_labels)==len(x) and bool(row.supported)==(len(x)>=3)
                expected=np.array([np.mean(x),*np.quantile(x,[.25,.75])]) if len(x) else np.full(3,np.nan)
                summary_error=max(summary_error,equal(row[['raw_mean','raw_q25','raw_q75']],expected))
                visible=expected if len(x)>=3 else np.full(3,np.nan)
                summary_error=max(summary_error,equal(row[['mean','q25','q75']],visible))
                support_count+=int(len(x)>=3)
        records.append(dict(path=stem,cells=len(assignment),library_error=library_error,
            specimen_error=specimen_error,reported_label_error=label_error,
            summary_error=summary_error,supported_feature_bins=support_count,
            within_reported_label_changes_checked=checked_changes,
            unsupported_patient_links_excluded_from_main_summary=True))
    assert {r['path'] for r in records}==expected_paths
    if len(methods) == 1:
        assert records
        return dict(paths=records,scope=methods[0],cross_method_comparison='pending other method',
                    all_library_thresholds_and_label_support_rules_verified=True)
    # Ensure that all method/path pairs were reported, including low overlap.
    comparison=pd.read_csv(TAB/'all_path_method_comparisons.csv')
    required=set()
    for root in ['primary','alternative']:
        a,b=source['Monocle2',root],source['Slingshot',root]
        for ac in [c for c in a if c.startswith('PT_')]:
            for bc in [c for c in b if c.startswith('PT_')]:
                required.add((root,ac,bc))
                valid=a[ac].notna()&b[bc].notna()&a['W_'+ac.split('_')[1]].ge(.5)&b['W_'+bc.split('_')[1]].ge(.5)
                row=comparison[comparison.root.eq(root)&comparison.monocle_path.eq(ac)&comparison.slingshot_path.eq(bc)]
                assert len(row)==1 and row.overlapping_supported_cells.iloc[0]==int(valid.sum())
                if valid.sum()>=30 and a.loc[valid,ac].nunique()>1 and b.loc[valid,bc].nunique()>1:
                    rank_a=a.loc[valid,ac].rank().to_numpy()
                    rank_b=b.loc[valid,bc].rank().to_numpy()
                    equal(row.descriptive_spearman,np.array([np.corrcoef(rank_a,rank_b)[0,1]]))
                else:
                    assert row.descriptive_spearman.isna().all()
    actual=set(comparison[['root','monocle_path','slingshot_path']].itertuples(index=False,name=None))
    assert actual==required and records
    return dict(paths=records,all_method_path_pairs_reported=True,pairs=len(required),
                all_library_thresholds_and_label_support_rules_verified=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['all','scores-and-compositions','dynamics'],default='all')
    parser.add_argument('--method',choices=['all','Monocle2','Slingshot'],default='all')
    args=parser.parse_args();stage=args.stage
    methods=['Monocle2','Slingshot'] if args.method=='all' else [args.method]
    metadata=pd.read_csv(ROOT/"work/AllEpithelial/metadata.csv.gz",low_memory=False)
    assert len(metadata)==223736 and metadata.cell_id.is_unique
    values=aligned(TAB/"function_and_CNV_values.csv.gz",metadata)
    if stage=='dynamics':
        report=json.loads((LOG/'function_and_composition_independent_validation.json').read_text())
        for path,digest in report['verified_output_hashes'].items():
            with (ROOT/path).open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==digest
        report['dynamics']=validate_dynamics(metadata,values,methods)
    else:
        report=dict(scores=validate_scores(metadata,values),compositions=validate_compositions(metadata),passed=True)
        files=['function_and_CNV_values.csv.gz','function_gene_standardization.csv','function_gene_coverage.csv','function_gene_membership.csv']
        report['verified_output_hashes']={}
        for name in files:
            path=TAB/name
            with path.open('rb') as f:report['verified_output_hashes'][str(path.relative_to(ROOT))]=hashlib.file_digest(f,'sha256').hexdigest()
        if stage=='all':report['dynamics']=validate_dynamics(metadata,values,methods)
    filename='function_and_composition_independent_validation.json' if stage=='scores-and-compositions' else 'postprocessing_independent_validation.json'
    if args.method!='all':filename=args.method+'_'+filename
    (LOG/filename).write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2),flush=True)


if __name__=="__main__":main()
