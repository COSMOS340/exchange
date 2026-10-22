"""Source composition and unit-equal descriptive profiles on the fixed tree."""
import json
import numpy as np
import pandas as pd

from pipeline_common import ROOT,BLOCK,CELLS,TAB,LOG,fixed_bins,equal_specimen_bins,save_json


def main():
    assert json.loads((LOG/'Monocle2_independent_validation.json').read_text())['passed']
    m=pd.read_csv(ROOT/'work'/BLOCK/'metadata.csv.gz',low_memory=False)
    assert len(m)==CELLS and m.cell_id.is_unique
    folder=TAB/'source_comparisons';folder.mkdir(exist_ok=False)
    traits={'tissue_group':['Normal_colorectal','Primary_CRC'],
            'display_3class':['Normal','Adjacent','CRC'],
            'CNV_relative_band':['Low','Intermediate','High'],
            'broad_identity':sorted(m.broad_identity.unique()),'dataset':sorted(m.dataset.unique())}
    features={field+'::'+str(level):(field,level) for field,levels in traits.items() for level in levels}
    rows=[]
    for identity,frame in m.groupby('broad_identity',sort=True):
        row=dict(broad_identity=identity,cells=len(frame))
        for source in traits['tissue_group']:
            d=frame[frame.tissue_group.eq(source)]
            row[source+'_cells']=len(d);row[source+'_captures']=d.capture_id.nunique()
        row['shared_display_eligible']=all(row[source+'_cells']>=30 and row[source+'_captures']>=3 for source in traits['tissue_group'])
        rows.append(row)
    pd.DataFrame(rows).to_csv(folder/'shared_identity_source_coverage.csv',index=False)
    manifest=[]
    for root in ['primary','alternative']:
        cells=pd.read_csv(TAB/('Monocle2_'+root+'_cells.csv.gz'))
        assert np.array_equal(cells.cell_id,m.cell_id)
        time=cells.PT_1.to_numpy();scaled=(time-time.min())/np.ptp(time)
        frame=m[['cell_id','dataset','subject_key','specimen_id','capture_id','tissue_group',
                 'display_3class','broad_identity','cnv_abs_log2','nCount_input']].copy()
        frame['bin']=fixed_bins(scaled,30)
        for feature,(field,level) in features.items():frame[feature]=m[field].eq(level).astype(float)
        captures,specimens=equal_specimen_bins(frame,list(features))
        captures.to_csv(folder/(root+'_source_capture_bins.csv.gz'),index=False)
        specimens.to_csv(folder/(root+'_source_specimen_bins.csv.gz'),index=False)
        annotations=[]
        for b in range(30):
            pooled=frame[frame.bin.eq(b)];cap=captures[captures.bin.eq(b)&captures.eligible]
            spec=specimens[specimens.bin.eq(b)]
            for feature,(field,level) in features.items():
                annotations.append(dict(bin=b,field=field,level=level,
                    pooled_cell_fraction=float(pooled[feature].mean()) if len(pooled) else np.nan,
                    equal_capture_fraction=float(cap[feature].mean()) if len(cap) else np.nan,
                    equal_capture_then_specimen_fraction=float(spec[feature].mean()) if len(spec) else np.nan,
                    cells=len(pooled),eligible_captures=len(cap),specimens=len(spec),datasets=spec.dataset.nunique()))
        pd.DataFrame(annotations).to_csv(folder/(root+'_dynamic_heatmap_annotations.csv.gz'),index=False)
        # Histograms are probability masses. Zero genuinely means no cells in a
        # bin; this is distinct from missing expression support.
        distribution=[]
        for source,d in frame.groupby('display_3class',sort=False):
            counts=d.groupby(['dataset','specimen_id','capture_id','bin']).size()
            cap=counts.unstack('bin',fill_value=0).reindex(columns=range(30),fill_value=0)
            cap=cap.div(cap.sum(axis=1),axis=0)
            spec=cap.groupby(level=[0,1]).mean()
            pooled=d.bin.value_counts().reindex(range(30),fill_value=0)/len(d)
            for b in range(30):distribution.append(dict(source=source,bin=b,pooled_mass=pooled[b],
                equal_capture_then_specimen_mass=spec[b].mean(),cells=len(d),captures=len(cap),specimens=len(spec)))
        pd.DataFrame(distribution).to_csv(folder/(root+'_source_PT_distributions.csv'),index=False)
        # Source-conditional and shared-identity CNV trends remain descriptive.
        profile=[]
        identities=['All_identities',*sorted(m.broad_identity.unique())]
        for identity in identities:
            d=frame if identity=='All_identities' else frame[frame.broad_identity.eq(identity)]
            for source,z in d.groupby('display_3class',sort=False):
                cap,spec=equal_specimen_bins(z,['cnv_abs_log2'])
                for b in range(30):
                    values=spec.loc[spec.bin.eq(b),'cnv_abs_log2'].to_numpy()
                    obs=z[z.bin.eq(b)]
                    profile.append(dict(root=root,identity=identity,source=source,bin=b,
                        pooled_cell_mean=obs.cnv_abs_log2.mean(),raw_specimen_equal_mean=values.mean() if len(values) else np.nan,
                        mean=values.mean() if len(values)>=3 else np.nan,
                        q25=np.quantile(values,.25) if len(values)>=3 else np.nan,
                        q75=np.quantile(values,.75) if len(values)>=3 else np.nan,
                        cells=len(obs),specimens=len(values),datasets=spec.loc[spec.bin.eq(b),'dataset'].nunique()))
        pd.DataFrame(profile).to_csv(folder/(root+'_source_CNV_profiles.csv.gz'),index=False)
        leaveout=[]
        for excluded in sorted(m.dataset.unique()):
            remain=specimens[specimens.dataset.ne(excluded)]
            for b in range(30):
                d=remain[remain.bin.eq(b)]
                for source in traits['display_3class']:
                    feature='display_3class::'+source
                    leaveout.append(dict(excluded_dataset=excluded,source=source,bin=b,
                        equal_specimen_fraction=d[feature].mean() if len(d)>=3 else np.nan,
                        specimens=len(d),datasets=d.dataset.nunique(),
                        graph_and_Harmony_and_curves_refitted=False))
        pd.DataFrame(leaveout).to_csv(folder/(root+'_leave_one_cohort_source_profiles.csv'),index=False)
        manifest.append(dict(root=root,cells=CELLS,unmatched_sources=0,
            shared_identity_display_minimum_cells_per_2class=30,
            shared_identity_display_minimum_captures_per_2class=3,
            source_CNV_minimum_specimens_per_bin=3,
            Adjacent_only_dataset='GSE245552',inference_tests_performed=False))
    save_json(LOG/'source_description_tables.json',dict(completed=True,roots=manifest,
        independent_validation_pending=True,interpretation='Frozen common graph; source composition and transverse states, not real-time progression'))
    print('Source composition, shared-identity and exhaustive cohort-leaveout tables saved.',flush=True)


if __name__=='__main__':main()
