"""Independently reconstruct every gene/unit mean directly from raw CSC data."""
import json

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from pipeline_common import ROOT,BLOCK,CELLS,TAB,OBJ,LOG,save_json


def accumulate_direct(destination,groups,weights,counts):
    use=groups>=0
    if not use.any():return
    unique,local=np.unique(groups[use],return_inverse=True)
    projection=sp.csc_matrix((weights[use],(np.flatnonzero(use),local)),
                             shape=(counts.shape[1],len(unique)))
    destination[unique]+=(counts@projection).T.toarray()


def main():
    m=pd.read_csv(ROOT/'work'/BLOCK/'metadata.csv.gz',low_memory=False)
    genes=pd.read_csv(ROOT/'work'/BLOCK/'genes.csv').gene.to_numpy()
    folder=TAB/'dynamic_genes'
    detection=pd.read_csv(folder/'dynamic_gene_detection_coverage.csv')
    assert np.array_equal(detection.gene,genes)
    layouts={}
    for root in ['primary','alternative']:
        cells=pd.read_csv(TAB/('Monocle2_'+root+'_cells.csv.gz'))
        assert np.array_equal(cells.cell_id,m.cell_id)
        pt=cells.PT_1.to_numpy();scaled=(pt-pt.min())/(pt.max()-pt.min())
        bins=np.minimum(np.searchsorted(np.linspace(0,1,31),scaled,side='right')-1,29)
        assignments=pd.read_csv(folder/('Monocle2_'+root+'_cell_bins.csv.gz'))
        assert np.array_equal(assignments.cell_id,m.cell_id)
        assert np.array_equal(assignments.bin,bins)
        capture=pd.read_csv(folder/('Monocle2_'+root+'_capture_bins.csv.gz'))
        specimen=pd.read_csv(folder/('Monocle2_'+root+'_specimen_bins.csv.gz'))
        assert not capture.duplicated(['capture_id','bin']).any()
        assert not specimen.duplicated(['specimen_id','bin']).any()
        cap_lookup={(row.capture_id,int(row.bin)):i for i,row in enumerate(capture.itertuples())}
        spec_lookup={(row.specimen_id,int(row.bin)):i for i,row in enumerate(specimen.itertuples())}
        cap_group=np.array([cap_lookup[(cap,int(b))] for cap,b in zip(m.capture_id,bins)])
        cap_count=np.bincount(cap_group,minlength=len(capture))
        assert np.array_equal(cap_count,capture.cells)
        assert np.array_equal(capture.eligible,capture.cells>=10)
        eligible=cap_count[cap_group]>=10
        spec_group=np.array([spec_lookup.get((spec,int(b)),-1) if good else -1
                             for spec,b,good in zip(m.specimen_id,bins,eligible)])
        n_captures=np.zeros(len(specimen),dtype=int)
        for row in capture[capture.eligible].itertuples():
            n_captures[spec_lookup[(row.specimen_id,int(row.bin))]]+=1
        assert np.array_equal(n_captures,specimen.captures)
        cap_weights=1/cap_count[cap_group]
        spec_weights=np.zeros(CELLS)
        spec_weights[eligible]=cap_weights[eligible]/n_captures[spec_group[eligible]]
        layouts[root]=dict(cap_group=cap_group,spec_group=spec_group,cap_weights=cap_weights,
            spec_weights=spec_weights,reconstructed_capture=np.zeros((len(capture),len(genes))),
            reconstructed_specimen=np.zeros((len(specimen),len(genes))))
    spec_codes,spec_ids=pd.factorize(m.specimen_id,sort=False)
    positive_by_specimen=np.zeros((len(spec_ids),len(genes)))
    log_sum,log_square=np.zeros(len(genes)),np.zeros(len(genes))
    with h5py.File(ROOT/'data/Overall/counts.h5') as h:
        ptr=h['indptr'][:]
        for lo in range(0,CELLS,3179):
            hi=min(lo+3179,CELLS);a,b=int(ptr[lo]),int(ptr[hi])
            raw=sp.csc_matrix((h['data'][a:b].astype(float),h['indices'][a:b],ptr[lo:hi+1]-a),
                               shape=(len(genes),hi-lo))
            assert np.array_equal(np.asarray(raw.sum(axis=0)).ravel(),m.common_gene_UMI.iloc[lo:hi])
            positive=raw.copy();positive.data=(positive.data>0).astype(float);positive.eliminate_zeros()
            accumulate_direct(positive_by_specimen,spec_codes[lo:hi],np.ones(hi-lo),positive)
            expression=raw.multiply((1e4/m.nCount_input.iloc[lo:hi].to_numpy())[None,:]).tocsc()
            expression.data=np.log1p(expression.data)
            log_sum+=np.asarray(expression.sum(axis=1)).ravel()
            log_square+=np.asarray(expression.power(2).sum(axis=1)).ravel()
            for layout in layouts.values():
                accumulate_direct(layout['reconstructed_capture'],layout['cap_group'][lo:hi],layout['cap_weights'][lo:hi],expression)
                accumulate_direct(layout['reconstructed_specimen'],layout['spec_group'][lo:hi],layout['spec_weights'][lo:hi],expression)
    assert np.array_equal(positive_by_specimen.sum(axis=0),detection.detected_cells)
    assert np.array_equal((positive_by_specimen>=5).sum(axis=0),detection.specimens_with_at_least_5_positive_cells)
    assert np.array_equal(detection.detection_eligible,
                          (detection.detected_cells/CELLS>=.01)&(detection.specimens_with_at_least_5_positive_cells>=10))
    checks={}
    standardization=pd.read_csv(folder/'all_gene_expression_standardization.csv')
    means=log_sum/CELLS
    sd=np.sqrt(np.maximum(0,log_square/CELLS-means**2))
    assert np.array_equal(standardization.gene,genes)
    assert np.allclose(standardization.mean_log1p_CP10K,means,atol=1e-9,rtol=1e-10)
    assert np.allclose(standardization.population_sd,sd,atol=1e-9,rtol=1e-10)
    for root,layout in layouts.items():
        with h5py.File(OBJ/('gene_expression_bins_'+root+'.h5')) as h:
            assert np.array_equal(h['gene'].asstr()[:],genes)
            error={}
            for unit in ['capture','specimen']:
                actual=h[unit+'_bin_means'][:]
                expected=layout['reconstructed_'+unit]
                assert actual.shape==expected.shape and np.isfinite(actual).all()
                maximum=float(np.abs(actual-expected).max())
                assert maximum<1e-9,(root,unit,maximum)
                error[unit+'_mean_max_error']=maximum
                error[unit+'_scalar_values_verified']=int(actual.size)
            checks[root]=error
    save_json(LOG/'dynamic_expression_independent_validation.json',dict(passed=True,cells=CELLS,
        all_measured_genes=len(genes),source='Original input CSC, direct cell weights, different I/O partition',
        all_capture_and_specimen_gene_means_verified=True,all_gene_detection_preconditions_verified=True,
        fixed_bin_assignments_verified=True,roots=checks))
    print(json.dumps(checks,indent=2),flush=True)


if __name__=='__main__':main()
