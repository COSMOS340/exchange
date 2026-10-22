"""Validate every real gene using an explicit full-dummy QR WLS design."""
import json,subprocess
import h5py
import numpy as np
import pandas as pd
from pipeline_common import ROOT,TAB,OBJ,LOG,ENV,save_json


def independent_basis(t):
    knots=[0,0,0,0,1/3,2/3,1,1,1,1]
    def b(i,degree):
        if degree==0:return ((t>=knots[i])&(t<knots[i+1])).astype(float)
        left=0 if knots[i+degree]==knots[i] else (t-knots[i])/(knots[i+degree]-knots[i])*b(i,degree-1)
        right=0 if knots[i+degree+1]==knots[i+1] else (knots[i+degree+1]-t)/(knots[i+degree+1]-knots[i+1])*b(i+1,degree-1)
        return left+right
    return np.column_stack([b(i,3) for i in range(1,6)])


def main():
    reports=[]
    for root in ['primary','alternative']:
        folder=TAB/('dynamic_genes/Monocle2_'+root)
        observed=pd.read_csv(TAB/('dynamic_genes/Monocle2_'+root+'_specimen_bins.csv.gz'))
        coverage=pd.read_csv(folder/'dynamic_bin_coverage.csv')
        eligible=observed.bin.isin(coverage.loc[coverage.supported,'bin'])
        supported=observed[eligible]
        counts=supported.groupby('specimen_id').size()
        use=eligible & observed.specimen_id.map(counts).ge(2)
        obs=observed[use].reset_index(drop=True)
        codes,units=pd.factorize(obs.specimen_id,sort=True)
        dummy=np.eye(len(units))[codes]
        basis=independent_basis((obs.bin.to_numpy()+.5)/30)
        design=np.column_stack([dummy,basis])
        weights=1/np.bincount(codes)[codes]
        weighted=design*np.sqrt(weights[:,None])
        q,r=np.linalg.qr(weighted,mode='reduced')
        assert np.linalg.matrix_rank(r)==design.shape[1]
        projector=np.linalg.solve(r,q.T)*np.sqrt(weights[None,:])
        clusters=sorted(obs.cluster_unit.unique());g=len(clusters)
        cluster_indices=[np.flatnonzero(obs.cluster_unit.eq(c).to_numpy()) for c in clusters]
        correction=g/(g-1)*(len(obs)-1)/(len(obs)-design.shape[1])
        truth=pd.read_csv(folder/'dynamic_gene_statistics.csv.gz')
        with h5py.File(OBJ/('gene_expression_bins_'+root+'.h5')) as h:
            y=h['specimen_bin_means'][:][use]
            assert np.array_equal(h['gene'].asstr()[:],truth.gene)
        statistic=np.full(y.shape[1],np.nan);beta=np.empty((5,y.shape[1]));ranks=np.zeros(y.shape[1],int)
        for lo in range(0,y.shape[1],137):
            hi=min(lo+137,y.shape[1]);values=y[:,lo:hi]
            coefficients=projector@values;beta[:,lo:hi]=coefficients[-5:]
            residual=values-design@coefficients
            cluster_influence=np.stack([projector[-5:,idx]@residual[idx] for idx in cluster_indices])
            covariance=np.einsum('gij,gkj->jik',cluster_influence,cluster_influence)*correction
            eigen=np.linalg.eigvalsh(covariance);rank=(eigen>np.maximum(np.maximum(eigen[:,-1],0)[:,None]*1e-10,1e-20)).sum(axis=1)
            ranks[lo:hi]=rank;assess=rank==5
            v=coefficients[-5:,assess].T
            statistic[np.arange(lo,hi)[assess]]=np.maximum(0,np.einsum('gi,gi->g',v,np.linalg.solve(covariance[assess],v[...,None])[...,0])/5)
        actual_beta=truth[[f'time_beta_{i}' for i in range(1,6)]].to_numpy().T
        assert np.allclose(beta,actual_beta,rtol=2e-7,atol=2e-9)
        assert np.array_equal(ranks,truth.covariance_rank)
        tested=truth.tested.to_numpy()
        assert np.allclose(statistic[tested],truth.F_statistic[tested],rtol=2e-7,atol=2e-8)
        numeric=LOG/f'{root}_dynamic_independent_F.csv'
        pd.DataFrame(dict(F_statistic=statistic,tested=tested)).to_csv(numeric,index=False)
        ref=LOG/f'{root}_dynamic_independent_R_probabilities.csv'
        subprocess.run([str(ENV/'bin/Rscript'),str(ROOT/'results/code/validate_dynamic_probabilities.R'),str(numeric),str(ref),str(g-1)],check=True)
        probabilities=pd.read_csv(ref)
        for field in ['p_value','BH_q_value']:
            assert np.allclose(probabilities[field],truth[field],rtol=2e-6,atol=2e-10,equal_nan=True),field
        curve=independent_basis((np.arange(30)+.5)/30)@beta
        amplitude=np.ptp(curve[coverage.supported],axis=0)
        assert np.allclose(amplitude,truth.adjusted_dynamic_amplitude,rtol=2e-7,atol=2e-9)
        expected=tested&probabilities.BH_q_value.le(.05).to_numpy()&(amplitude>=.25)
        assert np.array_equal(expected,truth.selected)
        with h5py.File(OBJ/f'dynamic_model_{root}.h5') as h:
            assert np.allclose(curve[coverage.supported],h['adjusted_curves'][:][coverage.supported],rtol=2e-7,atol=2e-9)
        reports.append(dict(root=root,all_genes=y.shape[1],all_tested=int(tested.sum()),
            full_dummy_columns=design.shape[1],specimens=len(units),observations=len(obs),clusters=g,
            maximum_beta_absolute_error=float(np.max(np.abs(beta-actual_beta))),
            maximum_F_absolute_error=float(np.max(np.abs(statistic[tested]-truth.F_statistic[tested]))),
            selected_identical=int(expected.sum()),explicit_QR_WLS_and_cluster_sandwich=True,R_pf_and_BH_checked=True))
        print('Validated actual full gene model',root,flush=True)
    save_json(LOG/'dynamic_model_full_independent_validation.json',dict(passed=True,roots=reports,cell_subsampling=0))


if __name__=='__main__':main()
