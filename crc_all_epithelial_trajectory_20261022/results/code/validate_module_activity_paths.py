"""Reconstruct every module cell score from original CSC and every path unit."""
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pipeline_common import ROOT,BLOCK,CELLS,TAB,LOG,save_json,digest


def same(a,b):
    assert np.allclose(np.asarray(a,float),np.asarray(b,float),rtol=2e-9,atol=2e-10,equal_nan=True)


def main():
    m=pd.read_csv(ROOT/'work'/BLOCK/'metadata.csv.gz',low_memory=False)
    folder=TAB/'dynamic_genes';source_path=folder/'module_and_representative_values_all_cells.csv.gz'
    scores=pd.read_csv(source_path);assert np.array_equal(scores.cell_id,m.cell_id)
    standard=pd.read_csv(folder/'all_gene_expression_standardization.csv')
    assignments=pd.read_csv(folder/'Monocle2_primary/dynamic_gene_module_assignment.csv')
    reps=pd.read_csv(folder/'Monocle2_primary/dynamic_module_representative_genes.csv')
    modules=sorted(assignments.loc[assignments.selected,'module'].unique())
    features=modules+reps.sort_values(['module','rank']).gene.tolist()
    indices={module:np.flatnonzero(assignments.module.eq(module).to_numpy()) for module in modules}
    by_gene={g:i for i,g in enumerate(standard.gene)}
    mean=standard.mean_log1p_CP10K.to_numpy();sd=standard.population_sd.to_numpy()
    errors={f:0. for f in features}
    with h5py.File(ROOT/'data/Overall/counts.h5') as h:
        ptr=h['indptr'][:];assert tuple(h['shape'][:])==(15635,CELLS)
        for lo in range(0,CELLS,3479):
            hi=min(CELLS,lo+3479);a,b=int(ptr[lo]),int(ptr[hi])
            matrix=sp.csc_matrix((h['data'][a:b],h['indices'][a:b],ptr[lo:hi+1]-a),shape=(15635,hi-lo))
            library=m.nCount_input.iloc[lo:hi].to_numpy()
            for module in modules:
                member=indices[module]
                dense=matrix[member,:].toarray().astype(float)
                log=np.log1p(dense*(1e4/library)[None,:])
                expected=((log-mean[member,None])/sd[member,None]).mean(axis=0)
                actual=scores[module].iloc[lo:hi].to_numpy();same(expected,actual)
                errors[module]=max(errors[module],float(np.max(np.abs(expected-actual))))
            for gene in features[len(modules):]:
                expected=np.log1p(matrix[[by_gene[gene]],:].toarray().ravel()*1e4/library)
                actual=scores[gene].iloc[lo:hi].to_numpy();same(expected,actual)
                errors[gene]=max(errors[gene],float(np.max(np.abs(expected-actual))))
    method_folder=folder/'method_path_sensitivity'
    summary=pd.read_csv(method_folder/'all_module_representative_path_summary.csv.gz')
    roster=pd.read_csv(method_folder/'all_method_root_path_roster.csv')
    records=[]
    for method in ['Monocle2','Slingshot']:
        for root in ['primary','alternative']:
            suffix='path_cells' if method=='Monocle2' else 'cells'
            paths=pd.read_csv(TAB/f'{method}_{root}_{suffix}.csv.gz');assert np.array_equal(paths.cell_id,m.cell_id)
            for column in [c for c in paths if c.startswith('PT_')]:
                number=column.split('_')[1];path='Path'+number
                pt=paths[column].to_numpy();valid=np.isfinite(pt)&paths['W_'+number].ge(.5).to_numpy()
                row=roster[roster.method.eq(method)&roster.root.eq(root)&roster.path.eq(path)]
                assert len(row)==1 and row.qualified_cells.iloc[0]==valid.sum()
                if row.status.iloc[0]!='descriptive':continue
                d=scores.loc[valid].copy();t=pt[valid]
                d['bin']=np.minimum(np.floor((t-t.min())/np.ptp(t)*12).astype(int),11)
                keys=['dataset','capture_id','specimen_id','subject_key','reported_patient_link_supported','bin']
                capture=d.groupby(keys)[features].mean();size=d.groupby(keys).size()
                qualified=capture[size.ge(10)]
                specimen=qualified.groupby(level=['dataset','specimen_id','subject_key','reported_patient_link_supported','bin']).mean()
                supported=specimen.index.get_level_values('reported_patient_link_supported')
                labels=specimen[supported].groupby(level=['dataset','subject_key','bin']).mean()
                proxy=specimen[~supported]
                actual=summary[summary.method.eq(method)&summary.root.eq(root)&summary.path.eq(path)]
                assert len(actual)==12*len(features)
                for record in actual.itertuples():
                    v=labels.loc[labels.index.get_level_values('bin')==record.bin,record.feature]
                    u=proxy.loc[proxy.index.get_level_values('bin')==record.bin,record.feature]
                    assert record.supported_reported_labels==len(v) and record.supported==(len(v)>=3)
                    same(record.raw_mean,v.mean() if len(v) else np.nan)
                    same([record.mean,record.q25,record.q75],[v.mean(),v.quantile(.25),v.quantile(.75)] if len(v)>=3 else [np.nan]*3)
                    assert record.unsupported_proxy_specimens==len(u)
                    same([record.proxy_mean,record.proxy_min,record.proxy_max],[u.mean(),u.min(),u.max()] if len(u) else [np.nan]*3)
                # Independently check every saved capture and specimen observation.
                name=f'{method}_{root}_{path}'
                saved=pd.read_csv(method_folder/(name+'_capture_bins.csv.gz')).set_index(keys).sort_index()
                same(saved[features],capture.sort_index());same(saved.cells,size.sort_index())
                saved_spec=pd.read_csv(method_folder/(name+'_specimen_bins.csv.gz')).set_index(['dataset','specimen_id','subject_key','reported_patient_link_supported','bin']).sort_index()
                same(saved_spec[features],specimen.sort_index())
                records.append(dict(method=method,root=root,path=path,cells=int(valid.sum()),capture_bins=len(capture),specimen_bins=len(specimen),summary_rows=len(actual)))
    assert len(roster)==sum(len([x for x in pd.read_csv(TAB/f'{method}_{root}_{"path_cells" if method=="Monocle2" else "cells"}.csv.gz',nrows=1) if x.startswith('PT_')]) for method in ['Monocle2','Slingshot'] for root in ['primary','alternative'])
    save_json(LOG/'module_activity_path_independent_validation.json',dict(passed=True,all_cells=CELLS,all_cell_features=len(features),
        direct_original_CSC_dense_member_gene_Z_max_errors=errors,all_returned_paths_checked=records,
        source_score_sha256=digest(source_path),all_path_masks_thresholds_aggregation_and_NA_checked=True))
    print('Every cell module/representative value and every returned path aggregation independently passed.',flush=True)


if __name__=='__main__':main()
