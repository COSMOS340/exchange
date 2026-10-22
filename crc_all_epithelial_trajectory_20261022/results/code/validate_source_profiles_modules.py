"""Independent source-profile arithmetic and full selected-module checks."""
import h5py
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from pipeline_common import ROOT,BLOCK,TAB,OBJ,LOG,save_json


def equal(a,b):
    a,b=np.asarray(a,float),np.asarray(b,float)
    assert np.allclose(a,b,rtol=2e-10,atol=2e-12,equal_nan=True),str((a,b))[:500]


def main():
    m=pd.read_csv(ROOT/'work'/BLOCK/'metadata.csv.gz',low_memory=False)
    folder=TAB/'source_comparisons';reports=[]
    for root in ['primary','alternative']:
        cells=pd.read_csv(TAB/f'Monocle2_{root}_cells.csv.gz');assert np.array_equal(cells.cell_id,m.cell_id)
        frame=m.copy();p=cells.PT_1.to_numpy();frame['bin']=np.minimum(np.floor((p-p.min())/np.ptp(p)*30).astype(int),29)
        annotations=pd.read_csv(folder/f'{root}_dynamic_heatmap_annotations.csv.gz')
        source_specs={}
        keys=['dataset','specimen_id','capture_id','bin']
        sizes=frame.groupby(keys).size()
        for field in annotations.field.unique():
            categorized=frame.assign(_category=frame[field])
            counts=categorized.groupby(keys+['_category'],observed=True).size().unstack('_category',fill_value=0)
            frac=counts.div(counts.sum(axis=1),axis=0);qualified=frac.loc[sizes.reindex(frac.index).ge(10)]
            spec=qualified.groupby(level=['dataset','specimen_id','bin']).mean()
            if field=='display_3class':source_specs=spec
            for row in annotations[annotations.field.eq(field)].itertuples():
                pooled=frame[frame.bin.eq(row.bin)]
                cap=qualified.loc[qualified.index.get_level_values('bin')==row.bin]
                s=spec.loc[spec.index.get_level_values('bin')==row.bin]
                equal(row.pooled_cell_fraction,pooled[field].eq(row.level).mean() if len(pooled) else np.nan)
                equal(row.equal_capture_fraction,cap[row.level].mean() if len(cap) else np.nan)
                equal(row.equal_capture_then_specimen_fraction,s[row.level].mean() if len(s) else np.nan)
                assert row.cells==len(pooled) and row.eligible_captures==len(cap) and row.specimens==len(s)
        leave=pd.read_csv(folder/f'{root}_leave_one_cohort_source_profiles.csv')
        for row in leave.itertuples():
            idx=source_specs.index
            subset=source_specs[(idx.get_level_values('dataset')!=row.excluded_dataset)&(idx.get_level_values('bin')==row.bin)]
            equal(row.equal_specimen_fraction,subset[row.source].mean() if len(subset)>=3 else np.nan)
            assert row.specimens==len(subset)
        distribution=pd.read_csv(folder/f'{root}_source_PT_distributions.csv')
        for source,d in frame.groupby('display_3class'):
            hist=d.groupby(['dataset','specimen_id','capture_id','bin']).size().unstack('bin',fill_value=0).reindex(columns=range(30),fill_value=0)
            cap_hist=hist.div(hist.sum(axis=1),axis=0)
            specimen_hist=cap_hist.groupby(level=['dataset','specimen_id']).mean()
            actual=distribution[distribution.source.eq(source)].sort_values('bin')
            equal(actual.pooled_mass,d.bin.value_counts(normalize=True).reindex(range(30),fill_value=0))
            equal(actual.equal_capture_then_specimen_mass,specimen_hist.mean(axis=0))
            equal(actual.equal_capture_then_specimen_mass.sum(),1)
        profiles=pd.read_csv(folder/f'{root}_source_CNV_profiles.csv.gz')
        for (identity,source),rows in profiles.groupby(['identity','source']):
            d=frame[frame.display_3class.eq(source)]
            if identity!='All_identities':d=d[d.broad_identity.eq(identity)]
            grouped=d.groupby(keys).cnv_abs_log2.agg(['mean','size'])
            values=grouped[grouped['size'].ge(10)]['mean'].groupby(level=['dataset','specimen_id','bin']).mean()
            for row in rows.itertuples():
                v=values[values.index.get_level_values('bin')==row.bin]
                obs=d[d.bin.eq(row.bin)]
                equal(row.pooled_cell_mean,obs.cnv_abs_log2.mean())
                equal(row.raw_specimen_equal_mean,v.mean() if len(v) else np.nan)
                equal([row['mean'] if isinstance(row,dict) else getattr(row,'mean'),row.q25,row.q75],
                    [v.mean(),v.quantile(.25),v.quantile(.75)] if len(v)>=3 else [np.nan]*3)
                assert row.specimens==len(v) and row.cells==len(obs)
        reports.append(dict(root=root,annotation_rows=len(annotations),leaveout_rows=len(leave),CNV_profile_rows=len(profiles),all_source_mass_rows=len(distribution)))
    old=ROOT.parent.parent/'exchange/crc_monocle2_20261021/work/CRC/metadata.csv.gz'
    literal=pd.read_csv(ROOT/'work'/BLOCK/'metadata.csv.gz',dtype=str,keep_default_na=False)
    literal=literal[literal.tissue_group.eq('Primary_CRC')].set_index('cell_id')
    previous=pd.read_csv(old,dtype=str,keep_default_na=False).set_index('cell_id').loc[literal.index]
    assert len(literal)==150917
    for field in ['CNV_relative_band','cnv_abs_log2','nCount_input','coarse_state','broad_identity','capture_id','specimen_id']:
        assert np.array_equal(literal[field],previous[field]),field
    a=pd.read_csv(TAB/'dynamic_genes/Monocle2_primary/dynamic_gene_module_assignment.csv');selected=a[a.selected]
    with h5py.File(OBJ/'dynamic_model_primary.h5') as h:
        z=h['gene_Z'][:].T[selected.index]
    use=np.isfinite(z).all(axis=0);z=z[:,use]
    dist=cdist(z,z);labels=selected.module.to_numpy();score=np.empty(len(z));centroids={}
    for module in sorted(set(labels)):
        mask=labels==module;within=dist[mask][:,mask].sum(axis=1)/(mask.sum()-1)
        between=np.min(np.stack([dist[mask][:,labels==other].mean(axis=1) for other in set(labels)-{module}]),axis=0)
        score[mask]=(between-within)/np.maximum(within,between)
        centroids[module]=z[mask].mean(axis=0)
    equal(score,selected.silhouette)
    metrics=pd.read_csv(TAB/'dynamic_genes/Monocle2_primary/dynamic_module_k_selection_metrics.csv')
    chosen=metrics[metrics.stable].sort_values(['silhouette','k'],ascending=[False,True]).iloc[0]
    assert chosen.k==len(set(labels))
    equal(chosen.silhouette,score.mean())
    reps=pd.read_csv(TAB/'dynamic_genes/Monocle2_primary/dynamic_module_representative_genes.csv')
    for module,d in reps.groupby('module'):
        mask=labels==module;genes=selected.loc[mask,'gene'].to_numpy();distance=np.square(z[mask]-centroids[module]).sum(axis=1)
        expected=genes[np.lexsort((genes,distance))[:3]]
        assert np.array_equal(expected,d.sort_values('rank').gene)
    save_json(LOG/'source_profiles_and_module_independent_validation.json',dict(passed=True,roots=reports,
        previous_CRC_150917_original_annotations_literal_identical=True,
        all_1787_gene_silhouettes_independently_recomputed=True,chosen_k=int(chosen.k),mean_silhouette=float(score.mean()),
        full_centroid_nearest_representatives_checked=True))
    print('Independent all source profiles, previous CRC labels and full module metrics passed.',flush=True)


if __name__=='__main__':main()
