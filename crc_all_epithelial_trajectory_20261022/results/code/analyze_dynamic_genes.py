"""Apply the frozen unit-based model; save every gene/bin including gaps."""
import json

import h5py
import numpy as np
import pandas as pd

from dynamic_model import bh_adjust, fit_dynamic_genes
from pipeline_common import ROOT, TAB, OBJ, LOG, save_json, validate_dynamic_source_tables


def bin_coverage(obs, bins):
    rows = []
    for b in range(bins):
        d = obs[obs.bin.eq(b)]
        rows.append(dict(bin=b, specimens=int(d.specimen_id.nunique()),
            captures=int(d.captures.sum()), cells=int(d.cells.sum()),
            datasets=int(d.dataset.nunique()),
            supported_reported_labels=int(d.loc[d.reported_patient_link_supported, 'subject_key'].nunique()),
            unsupported_proxy_specimens=int(d.loc[~d.reported_patient_link_supported, 'specimen_id'].nunique())))
    result = pd.DataFrame(rows)
    result['supported'] = result.specimens.ge(5) & result.datasets.ge(2)
    return result


def statistics_for(obs, y, coverage, genes, detection, bins):
    supported = coverage.supported.to_numpy()
    keep = obs.bin.map(coverage.set_index('bin').supported).to_numpy()
    model = fit_dynamic_genes(obs.loc[keep].reset_index(drop=True), y[keep], bins=bins)
    curves = model['adjusted_curves']
    amplitude = np.ptp(curves[supported], axis=0)
    tested = detection.detection_eligible.to_numpy() & np.isfinite(model['p'])
    raw_p = np.where(tested, model['p'], np.nan)
    q = bh_adjust(raw_p)
    support_ids = np.flatnonzero(supported)
    delta = curves[support_ids[-3:]].mean(axis=0)-curves[support_ids[:3]].mean(axis=0)
    table = detection.copy()
    table['tested'] = tested
    table['covariance_rank'] = model['covariance_rank']
    table['F_statistic'] = np.where(tested, model['statistic'], np.nan)
    table['p_value'] = raw_p
    table['BH_q_value'] = q
    table['adjusted_dynamic_amplitude'] = amplitude
    table['late_minus_early_adjusted_effect'] = delta
    table['direction'] = np.where(np.abs(delta) < .1*np.maximum(amplitude, 1e-12),
                                   'nonmonotonic_or_small_net_change', np.where(delta>0, 'net_increase', 'net_decrease'))
    table['selected'] = tested & (q<=.05) & (amplitude>=.25)
    for j in range(5):
        table['time_beta_' + str(j+1)] = model['beta'][j]
    model['curves_with_support'] = np.where(supported[:, None], curves, np.nan)
    return table, model


def analyze_root(root, bins):
    folder = TAB / 'dynamic_genes'
    output = folder / ('Monocle2_' + root)
    output.mkdir(exist_ok=False)
    obs = pd.read_csv(folder / ('Monocle2_' + root + '_specimen_bins.csv.gz'))
    detection = pd.read_csv(folder / 'dynamic_gene_detection_coverage.csv')
    with h5py.File(OBJ / ('gene_expression_bins_' + root + '.h5')) as h:
        genes = h['gene'].asstr()[:]
        y = h['specimen_bin_means'][:]
    assert np.array_equal(genes, detection.gene)
    coverage = bin_coverage(obs, bins)
    if coverage.supported.sum()<6:
        raise RuntimeError('Fewer than six supported bins; record unassessable dynamic model instead of changing thresholds')
    coverage.to_csv(output / 'dynamic_bin_coverage.csv', index=False)
    statistics, model = statistics_for(obs, y, coverage, genes, detection, bins)
    mean = np.full((bins, len(genes)), np.nan)
    q25, q75 = mean.copy(), mean.copy()
    for b in range(bins):
        values = y[obs.bin.eq(b).to_numpy()]
        if len(values):
            mean[b] = values.mean(axis=0)
            q25[b], q75[b] = np.quantile(values,[.25,.75],axis=0)
    curves = model['curves_with_support']
    centered = curves-np.nanmean(curves,axis=0)
    sd = np.nanstd(curves,axis=0)
    z = np.divide(centered, sd[None,:], out=np.full_like(centered,np.nan), where=sd[None,:]>1e-12)
    source_keys = dict(gene=np.repeat(genes,bins),bin=np.tile(np.arange(bins),len(genes)))
    binned = pd.DataFrame(source_keys)
    binned['raw_specimen_equal_mean_log1p_CP10K'] = mean.T.ravel()
    binned['raw_specimen_q25'] = q25.T.ravel()
    binned['raw_specimen_q75'] = q75.T.ravel()
    binned['supported'] = np.tile(coverage.supported.to_numpy(),len(genes))
    binned['display_mean'] = np.where(binned.supported,binned.raw_specimen_equal_mean_log1p_CP10K,np.nan)
    smooth = pd.DataFrame(source_keys)
    smooth['adjusted_spline_effect_log1p_CP10K'] = curves.T.ravel()
    smooth['gene_Z'] = z.T.ravel()
    smooth['supported'] = binned.supported
    validate_dynamic_source_tables(statistics,binned,smooth,bins)
    statistics.to_csv(output / 'dynamic_gene_statistics.csv.gz',index=False)
    binned.to_csv(output / 'dynamic_gene_binned_expression.csv.gz',index=False)
    smooth.to_csv(output / 'dynamic_gene_smoothed_expression.csv.gz',index=False)
    sensitivity=[]
    for q in [.01,.05,.1]:
        for amplitude in [.1,.25,.5]:
            selected=statistics.tested & statistics.BH_q_value.le(q) & statistics.adjusted_dynamic_amplitude.ge(amplitude)
            sensitivity.append(dict(q=q,amplitude=amplitude,genes=int(selected.sum())))
    pd.DataFrame(sensitivity).to_csv(output/'dynamic_threshold_sensitivity.csv',index=False)
    with h5py.File(OBJ/('dynamic_model_' + root + '.h5'),'x') as h:
        h.create_dataset('gene',data=genes.astype(object),dtype=h5py.string_dtype())
        h.attrs['expression_source'] = 'gene_expression_bins_' + root + '.h5'
        h.create_dataset('beta',data=model['beta'])
        h.create_dataset('adjusted_curves',data=curves)
        h.create_dataset('gene_Z',data=z)
    summary={k:v for k,v in model.items() if k not in ['beta','adjusted_curves','curves_with_support','statistic','p','covariance_rank']}
    summary.update(root=root,measured_genes=len(genes),detection_eligible=int(detection.detection_eligible.sum()),
                   tested_genes=int(statistics.tested.sum()),selected_genes=int(statistics.selected.sum()),
                   supported_bins=int(coverage.supported.sum()),bins=bins,
                   source_tables_complete=True,independent_numeric_validation_pending=True,
                   inference_uses_cells_as_replicates=False)
    save_json(output/'dynamic_model_summary.json',summary)
    print(root,json.dumps(summary),flush=True)


def main():
    bins=json.loads((ROOT/'ANALYSIS_PLAN.json').read_text())['dynamic_genes']['bins']
    if not json.loads((LOG/'dynamic_model_software_tests.json').read_text())['passed']:
        raise RuntimeError('Independent software tests have not passed')
    for root in ['primary','alternative']:
        analyze_root(root,bins)


if __name__=='__main__':
    main()
