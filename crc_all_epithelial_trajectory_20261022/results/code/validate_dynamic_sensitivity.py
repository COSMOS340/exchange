"""Validate every leaveout time coefficient with a separate full dummy QR design."""
import json
import h5py
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from pipeline_common import TAB,OBJ,LOG,save_json
from validate_dynamic_model_full import independent_basis


def ari(a,b):
    table=pd.crosstab(a,b).to_numpy(dtype=float)
    choose=lambda x:x*(x-1)/2
    n=table.sum();total=choose(n)
    observed=choose(table).sum();rows=choose(table.sum(axis=1)).sum();cols=choose(table.sum(axis=0)).sum()
    expected=rows*cols/total;denominator=(rows+cols)/2-expected
    return float((observed-expected)/denominator) if denominator else 1.


def main():
    folder=TAB/'dynamic_genes';sensitivity=folder/'sensitivity'
    roster=pd.read_csv(sensitivity/'all_cohort_leaveout_summary.csv')
    obs=pd.read_csv(folder/'Monocle2_primary_specimen_bins.csv.gz')
    with h5py.File(OBJ/'gene_expression_bins_primary.h5') as h:y=h['specimen_bin_means'][:]
    baseline=pd.read_csv(folder/'Monocle2_primary/dynamic_gene_module_assignment.csv')
    original=set(baseline.loc[baseline.selected,'gene'])
    reports=[]
    for record in roster.itertuples():
        if record.status!='assessed':
            assert (sensitivity/('leaveout_'+record.excluded_dataset)/'unassessable.json').exists();continue
        dest=sensitivity/('leaveout_'+record.excluded_dataset)
        cov=pd.read_csv(dest/'bin_coverage.csv')
        use=~obs.dataset.eq(record.excluded_dataset)&obs.bin.isin(cov.loc[cov.supported,'bin'])
        counts=obs[use].groupby('specimen_id').size();use&=obs.specimen_id.map(counts).ge(2)
        selected=obs[use].reset_index(drop=True)
        code,units=pd.factorize(selected.specimen_id,sort=True)
        x=np.column_stack([np.eye(len(units))[code],independent_basis((selected.bin.to_numpy()+.5)/30)])
        weight=1/np.bincount(code)[code];q,r=np.linalg.qr(x*np.sqrt(weight[:,None]),mode='reduced')
        projection=np.linalg.solve(r,q.T)[-5:]*np.sqrt(weight[None,:])
        expected=projection@y[use]
        table=pd.read_csv(dest/'all_gene_statistics.csv.gz')
        actual=table[[f'time_beta_{j}' for j in range(1,6)]].to_numpy().T
        assert np.allclose(actual,expected,rtol=2e-7,atol=2e-9)
        curve=independent_basis((np.arange(30)+.5)/30)@expected
        support=cov.supported.to_numpy();amplitude=np.ptp(curve[support],axis=0)
        assert np.allclose(amplitude,table.adjusted_dynamic_amplitude,rtol=2e-7,atol=2e-9)
        values=set(table.loc[table.selected,'gene'])
        assert record.selected_genes==len(values) and record.intersection_primary==len(original&values)
        assert abs(record.selection_Jaccard-len(original&values)/len(original|values))<1e-12
        partition=pd.read_csv(dest/'primary_gene_partition.csv')
        assert np.array_equal(partition.gene,baseline.loc[baseline.selected,'gene'])
        expected_ari=ari(partition.fixed_primary_module,partition.leaveout_fixed_k_label)
        assert abs(expected_ari-record.primary_genes_partition_ARI)<1e-12
        normalized=curve[support][:,baseline.selected].T
        normalized=(normalized-normalized.mean(axis=1)[:,None])/normalized.std(axis=1)[:,None]
        # Every final member is assigned to a nearest fitted centroid in leaveout profile space.
        labels=partition.leaveout_fixed_k_label.to_numpy();centers=np.array([normalized[labels==c].mean(axis=0) for c in range(labels.max()+1)])
        distance=cdist(normalized,centers,metric='sqeuclidean')
        assert np.all(distance[np.arange(len(labels)),labels]<=distance.min(axis=1)+1e-7)
        reports.append(dict(excluded_dataset=record.excluded_dataset,genes=len(table),specimens=len(units),
            all_five_time_coefficients_max_error=float(np.max(np.abs(expected-actual))),all_amplitudes_and_partition_ARI_checked=True))
        print('Independent leaveout coefficients',record.excluded_dataset,flush=True)
    root=pd.read_csv(sensitivity/'root_gene_selection_comparison.csv.gz')
    report=json.loads((LOG/'dynamic_sensitivity_summary.json').read_text())
    p=set(root.loc[root.selected_primary,'gene']);a=set(root.loc[root.selected_alternative,'gene'])
    assert len(p&a)==report['root_selection_overlap']
    assert abs(len(p&a)/len(p|a)-report['root_selection_Jaccard'])<1e-12
    save_json(LOG/'dynamic_sensitivity_independent_validation.json',dict(passed=True,
        every_cohort_attempt_checked=len(roster)==8,assessed=reports,root_selection_overlap_checked=True,
        inference_model_software_and_primary_alternative_full_F_P_BH_checks='dynamic_model_full_independent_validation.json',
        leaveout_F_P_covariance_method_reuses_the_independently_validated_model=True))


if __name__=='__main__':main()
