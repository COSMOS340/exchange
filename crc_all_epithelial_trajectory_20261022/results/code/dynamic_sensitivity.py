"""Fixed-graph root and all-cohort leaveout sensitivities; no favorable selection."""
import json
import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from analyze_dynamic_genes import bin_coverage,statistics_for
from cluster_profiles import standardize_profiles,cluster_all_profiles,adjusted_rand
from pipeline_common import ROOT,BLOCK,TAB,OBJ,LOG,save_json


def roster():
    m=pd.read_csv(ROOT/'work'/BLOCK/'metadata.csv.gz',low_memory=False)
    base=m.groupby('specimen_id',sort=True).agg(dataset=('dataset','first'),cells=('cell_id','size'),
        captures=('capture_id','nunique'),subject_key=('subject_key','first'),
        reported_patient_link_supported=('reported_patient_link_supported','first'),
        tissue_group=('tissue_group','first'),display_3class=('display_3class','first'))
    assert len(base)==229 and int(base.cells.sum())==223736
    reports=[]
    for root in ['primary','alternative']:
        captures=pd.read_csv(TAB/f'dynamic_genes/Monocle2_{root}_capture_bins.csv.gz')
        obs=pd.read_csv(TAB/f'dynamic_genes/Monocle2_{root}_specimen_bins.csv.gz')
        cover=pd.read_csv(TAB/f'dynamic_genes/Monocle2_{root}/dynamic_bin_coverage.csv')
        inf=obs[obs.bin.isin(cover.loc[cover.supported,'bin'])]
        d=base.copy()
        d['observed_capture_bins']=captures.groupby('specimen_id').size().reindex(d.index,fill_value=0)
        d['eligible_capture_bins']=captures[captures.eligible].groupby('specimen_id').size().reindex(d.index,fill_value=0)
        d['eligible_specimen_bins']=obs.groupby('specimen_id').size().reindex(d.index,fill_value=0)
        d['supported_specimen_bins']=inf.groupby('specimen_id').size().reindex(d.index,fill_value=0)
        d['model_informative']=d.supported_specimen_bins.ge(2)
        d['inference_status']=np.select([d.eligible_capture_bins.eq(0),d.supported_specimen_bins.eq(0),d.supported_specimen_bins.eq(1)],
            ['no_capture_bin_reaches_10_cells','only_globally_unsupported_bins','singleton_supported_bin'],default='informative_within_specimen')
        d['all_cells_retained_in_common_tree']=True
        d.to_csv(TAB/f'dynamic_genes/Monocle2_{root}_all_specimen_coverage.csv')
        reports.append(dict(root=root,total_specimens=len(d),total_cells=int(d.cells.sum()),
            status_counts=d.inference_status.value_counts().to_dict(),model_informative=int(d.model_informative.sum())))
    save_json(LOG/'all_specimen_dynamic_coverage.json',dict(passed=True,roots=reports,
        qualification_changes_only_downstream_summary_or_inference_not_common_tree=True))


def main():
    roster()
    folder=TAB/'dynamic_genes';out=folder/'sensitivity';out.mkdir(exist_ok=False)
    a=pd.read_csv(folder/'Monocle2_primary/dynamic_gene_module_assignment.csv')
    primary=pd.read_csv(folder/'Monocle2_primary/dynamic_gene_statistics.csv.gz')
    alternative=pd.read_csv(folder/'Monocle2_alternative/dynamic_gene_statistics.csv.gz')
    assert np.array_equal(primary.gene,alternative.gene)
    k=json.loads((folder/'Monocle2_primary/dynamic_module_selection.json').read_text())['chosen_k']
    comparison=primary[['gene','selected','BH_q_value','adjusted_dynamic_amplitude']].merge(
        alternative[['gene','selected','BH_q_value','adjusted_dynamic_amplitude']],on='gene',suffixes=('_primary','_alternative'),validate='one_to_one')
    comparison.to_csv(out/'root_gene_selection_comparison.csv.gz',index=False)
    pset=set(primary.loc[primary.selected,'gene']);aset=set(alternative.loc[alternative.selected,'gene'])
    rootstats=[]
    for root in ['primary','alternative']:
        cov=pd.read_csv(folder/f'Monocle2_{root}/dynamic_bin_coverage.csv')
        with h5py.File(OBJ/f'dynamic_model_{root}.h5') as h:curve=h['gene_Z'][:]
        rows=[]
        for module,d in a[a.selected].groupby('module'):
            profile=np.nanmean(curve[:,d.index],axis=1)
            rows.extend(dict(root=root,module=module,bin=b,mean_gene_Z=float(v) if np.isfinite(v) else np.nan,
                             primary_module_genes=len(d)) for b,v in enumerate(profile))
        pd.DataFrame(rows).to_csv(out/f'{root}_fixed_primary_module_profiles.csv',index=False)
        if root=='alternative':
            selected=alternative.selected.to_numpy()
            fit=cluster_all_profiles(standardize_profiles(curve[cov.supported][:,selected].T),k,20261023)
            cross=a.loc[selected,['gene','module','selected']].copy();cross['alternative_fixed_k_label']=fit['labels']
            cross.to_csv(out/'alternative_fixed_k_assignment.csv',index=False)
            common=cross.selected
            # a.selected is the primary-selection flag, including on alternate-selected rows.
            rootstats.append(dict(fixed_k=k,common_selected=int(common.sum()),
                ARI_common_genes=adjusted_rand(pd.factorize(cross.loc[common,'module'])[0],cross.loc[common,'alternative_fixed_k_label'].to_numpy())))
    p=pd.read_csv(TAB/'Monocle2_primary_cells.csv.gz');alt=pd.read_csv(TAB/'Monocle2_alternative_cells.csv.gz')
    assert np.array_equal(p.cell_id,alt.cell_id)
    orientation=float(spearmanr(p.PT_1,alt.PT_1).statistic)
    detection=pd.read_csv(folder/'dynamic_gene_detection_coverage.csv')
    obs=pd.read_csv(folder/'Monocle2_primary_specimen_bins.csv.gz')
    with h5py.File(OBJ/'gene_expression_bins_primary.h5') as h:y=h['specimen_bin_means'][:];genes=h['gene'].asstr()[:]
    baseline_labels=pd.factorize(a.loc[a.selected,'module'])[0]
    loo=[];profiles=[]
    for index,dataset in enumerate(sorted(obs.dataset.unique())):
        dest=out/('leaveout_'+dataset);dest.mkdir()
        try:
            use=~obs.dataset.eq(dataset).to_numpy();subset=obs.loc[use].reset_index(drop=True)
            coverage=bin_coverage(subset,30)
            if coverage.supported.sum()<6:raise ValueError('Fewer than six supported bins after leaveout')
            stats,model=statistics_for(subset,y[use],coverage,genes,detection,30)
            stats.to_csv(dest/'all_gene_statistics.csv.gz',index=False);coverage.to_csv(dest/'bin_coverage.csv',index=False)
            supported=coverage.supported.to_numpy();curves=model['curves_with_support']
            z=standardize_profiles(curves[supported][:,a.selected].T)
            fit=cluster_all_profiles(z,k,20261024+index)
            pd.DataFrame(dict(gene=a.loc[a.selected,'gene'],fixed_primary_module=a.loc[a.selected,'module'],
                leaveout_fixed_k_label=fit['labels'])).to_csv(dest/'primary_gene_partition.csv',index=False)
            for module,d in a[a.selected].groupby('module'):
                standardized=standardize_profiles(curves[supported][:,d.index].T).mean(axis=0)
                for b,value in zip(np.flatnonzero(supported),standardized):profiles.append(dict(excluded_dataset=dataset,module=module,bin=int(b),mean_gene_Z=float(value)))
            selected=set(stats.loc[stats.selected,'gene'])
            loo.append(dict(excluded_dataset=dataset,status='assessed',observations=model['observations'],specimens=model['specimens'],clusters=model['linkage_clusters'],
                supported_bins=int(supported.sum()),selected_genes=len(selected),intersection_primary=len(pset&selected),
                selection_Jaccard=len(pset&selected)/len(pset|selected),primary_genes_partition_ARI=adjusted_rand(baseline_labels,fit['labels'])))
        except (ValueError,np.linalg.LinAlgError) as e:
            save_json(dest/'unassessable.json',dict(excluded_dataset=dataset,status='unassessable',reason=str(e)))
            loo.append(dict(excluded_dataset=dataset,status='unassessable',reason=str(e)))
        print('Dynamic fixed-graph cohort sensitivity',dataset,loo[-1]['status'],flush=True)
    pd.DataFrame(loo).to_csv(out/'all_cohort_leaveout_summary.csv',index=False)
    pd.DataFrame(profiles).to_csv(out/'all_cohort_leaveout_fixed_module_profiles.csv',index=False)
    save_json(LOG/'dynamic_sensitivity_summary.json',dict(completed=True,primary_selected=len(pset),alternative_selected=len(aset),
        root_selection_overlap=len(pset&aset),root_selection_Jaccard=len(pset&aset)/len(pset|aset),
        root_cell_id_pseudotime_Spearman=orientation,alternative_partition=rootstats,
        all_eight_cohorts_attempted=len(loo)==8,cohort_results=loo,
        no_graph_Harmony_or_curves_refitted=True,independent_numeric_validation_pending=True))


if __name__=='__main__':main()
