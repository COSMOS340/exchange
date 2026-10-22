"""Choose gene modules using all-gene silhouette and whole-unit bootstrap."""
import json

import h5py
import numpy as np
import pandas as pd

from cluster_profiles import standardize_profiles,cluster_all_profiles,exact_silhouettes,adjusted_rand
from dynamic_model import prepare_design,spline_basis
from pipeline_common import ROOT,TAB,OBJ,LOG,save_json


def bootstrap_profiles(obs,y,bins,support,rng):
    """Resample linkage clusters within cohorts, not fitted biological cells."""
    groups=obs[['dataset','cluster_unit']].drop_duplicates()
    if groups.groupby('cluster_unit').dataset.nunique().gt(1).any():
        raise ValueError('Linkage cluster spans datasets unexpectedly')
    pieces=[];indices=[];replicate=0
    for dataset,stratum in groups.groupby('dataset',sort=True):
        labels=stratum.cluster_unit.to_numpy()
        for selected in rng.choice(labels,size=len(labels),replace=True):
            idx=np.flatnonzero(obs.cluster_unit.eq(selected).to_numpy())
            piece=obs.iloc[idx].copy()
            suffix='::bootstrap_copy_'+str(replicate)
            piece['specimen_id']=piece.specimen_id+suffix
            piece['cluster_unit']=piece.cluster_unit+suffix
            pieces.append(piece);indices.extend(idx);replicate+=1
    boot=pd.concat(pieces,ignore_index=True)
    values=y[np.asarray(indices)]
    d=prepare_design(boot,bins=bins)
    values=values[d['keep']]
    means=np.zeros((len(d['units']),values.shape[1]))
    np.add.at(means,d['unit'],values)
    means/=d['unit_counts'][:,None]
    centered=values-means[d['unit']]
    beta=d['bread']@(d['x'].T@(d['weights'][:,None]*centered))
    curve=spline_basis(d['bin_grid'])@beta
    return standardize_profiles(curve[support].T)


def main():
    folder=TAB/'dynamic_genes/Monocle2_primary'
    if (folder/'dynamic_module_selection.json').exists():
        raise RuntimeError('Module selection exists; preserve/version before recomputing')
    if not json.loads((LOG/'module_clustering_software_tests.json').read_text())['passed']:
        raise RuntimeError('Independent clustering metric tests have not passed')
    stats=pd.read_csv(folder/'dynamic_gene_statistics.csv.gz')
    selected=stats.selected.to_numpy()
    if selected.sum()<3:
        raise RuntimeError('Insufficient selected genes for modules; report without manufacturing a heatmap')
    coverage=pd.read_csv(folder/'dynamic_bin_coverage.csv')
    support=coverage.supported.to_numpy()
    with h5py.File(OBJ/'dynamic_model_primary.h5') as h:
        assert np.array_equal(h['gene'].asstr()[:],stats.gene)
        profiles=h['gene_Z'][:][support][:,selected].T
    if not np.isfinite(profiles).all():raise ValueError('Selected gene profiles are not finite')
    with h5py.File(OBJ/'gene_expression_bins_primary.h5') as h:
        y=h['specimen_bin_means'][:,selected]
    obs=pd.read_csv(TAB/'dynamic_genes/Monocle2_primary_specimen_bins.csv.gz')
    keep=obs.bin.map(coverage.set_index('bin').supported).to_numpy()
    obs=obs.loc[keep].reset_index(drop=True);y=y[keep]
    candidates={};skipped=[]
    for k in range(2,9):
        if k>=len(profiles):
            skipped.append(dict(k=k,reason='Fewer genes than requested modules'));continue
        try:candidates[k]=cluster_all_profiles(profiles,k,20261022+k*1000)
        except ValueError as e:skipped.append(dict(k=k,reason=str(e)))
    if not candidates:raise RuntimeError('No assessable module count')
    silhouettes=exact_silhouettes(profiles,{k:v['labels'] for k,v in candidates.items()})
    stability=[]
    for iteration in range(10):
        try:
            boot=bootstrap_profiles(obs,y,30,support,np.random.default_rng(20261022+iteration))
            for k,baseline in candidates.items():
                try:
                    fitted=cluster_all_profiles(boot,k,20271022+iteration*1000+k*10)
                    stability.append(dict(iteration=iteration,k=k,adjusted_Rand=adjusted_rand(baseline['labels'],fitted['labels']),status='assessed'))
                except ValueError as e:stability.append(dict(iteration=iteration,k=k,adjusted_Rand=np.nan,status=str(e)))
        except ValueError as e:
            for k in candidates:stability.append(dict(iteration=iteration,k=k,adjusted_Rand=np.nan,status=str(e)))
        print('Module unit-bootstrap iteration',iteration+1,'/ 10 complete',flush=True)
    stability=pd.DataFrame(stability)
    stability.to_csv(folder/'dynamic_module_bootstrap_stability.csv',index=False)
    metrics=[]
    for k,v in candidates.items():
        values=stability.loc[stability.k.eq(k),'adjusted_Rand'].dropna()
        median=float(values.median()) if len(values) else -1.
        metrics.append(dict(k=k,inertia=v['inertia'],silhouette=float(silhouettes[k].mean()),
                            median_bootstrap_ARI=median,valid_bootstrap_runs=len(values),
                            stable=bool(len(values)==10 and median>=.75)))
    metrics=pd.DataFrame(metrics)
    metrics.to_csv(folder/'dynamic_module_k_selection_metrics.csv',index=False)
    stable=metrics[metrics.stable]
    if len(stable):
        chosen=stable.sort_values(['silhouette','k'],ascending=[False,True]).iloc[0]
        decision='Maximum silhouette among stable candidates'
    else:
        chosen=metrics.sort_values(['median_bootstrap_ARI','silhouette','k'],ascending=[False,False,True]).iloc[0]
        decision='No candidate met stability criterion; exploratory partition chosen by stability then silhouette'
    k=int(chosen.k);fit=candidates[k]
    label_order=sorted(range(k),key=lambda c:(int(np.argmax(fit['centers'][c])),c))
    names={c:'Module'+str(i+1) for i,c in enumerate(label_order)}
    gene_names=stats.loc[selected,'gene'].to_numpy()
    assignments=stats[['gene','selected','BH_q_value','adjusted_dynamic_amplitude']].copy()
    assignments['module']='Not_selected'
    assignments.loc[selected,'module']=[names[c] for c in fit['labels']]
    assignments['silhouette']=np.nan
    assignments.loc[selected,'silhouette']=silhouettes[k]
    representatives=[];centers=[]
    for c in label_order:
        rows=np.flatnonzero(fit['labels']==c)
        distances=np.square(profiles[rows]-fit['centers'][c]).sum(axis=1)
        order=np.lexsort((gene_names[rows],distances))
        for rank,index in enumerate(rows[order[:3]]):
            representatives.append(dict(module=names[c],gene=gene_names[index],rank=rank+1,
                                        squared_centroid_distance=float(np.square(profiles[index]-fit['centers'][c]).sum())))
        centers.extend(dict(module=names[c],bin=int(b),mean_gene_Z=float(z),genes=len(rows))
                       for b,z in zip(np.flatnonzero(support),fit['centers'][c]))
    assignments.to_csv(folder/'dynamic_gene_module_assignment.csv',index=False)
    pd.DataFrame(representatives).to_csv(folder/'dynamic_module_representative_genes.csv',index=False)
    pd.DataFrame(centers).to_csv(folder/'dynamic_module_centroid_profiles.csv',index=False)
    save_json(folder/'dynamic_module_selection.json',dict(
        completed=True,chosen_k=k,selected_genes=int(selected.sum()),decision=decision,
        stability_threshold_met=bool(chosen.stable),bootstrap_runs=10,
        kmeans_initializations=5,kmeans_iterations_per_initialization=100,
        exact_silhouette_all_selected_genes=True,random_cells_sampled=0,
        bootstrap_unit='Whole reported linkage clusters, unsupported links pooled per dataset; stratified by dataset',
        selection_used_GO_or_tissue_CNV=False,skipped_candidates=skipped,
        independent_validation_pending=True))
    print('Selected',k,'modules;',decision,flush=True)


if __name__=='__main__':main()
