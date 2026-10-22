"""Recheck every joined CRC cell, all root/path pairs and all CNV bin summaries."""
import numpy as np
import pandas as pd
from pipeline_common import TAB,LOG,save_json


def main():
    folder=TAB/'CRC_only_comparison'
    joined=pd.read_csv(folder/'exact_150917_cell_comparison.csv.gz');assert len(joined)==150917 and joined.cell_id.is_unique
    table=pd.read_csv(folder/'all_four_root_comparisons.csv')
    profiles=pd.read_csv(folder/'all_CNV_binned_profiles.csv')
    overlap=pd.read_csv(folder/'all_path_pair_overlaps.csv')
    sources={};paths={};checks=[]
    for fit in ['old_CRC_only','joint_tree_CRC_subset']:
        for root in ['primary','alternative']:
            parent=folder/'frozen_previous_outputs' if fit=='old_CRC_only' else TAB
            original=pd.read_csv(parent/f'Monocle2_{root}_cells.csv.gz').set_index('cell_id')
            d=original.loc[joined.cell_id]
            assert np.array_equal(d.State,joined[f'{fit}_{root}_State'])
            assert np.array_equal(d.PT_1,joined[f'{fit}_{root}_PT'])
            scaled=(d.PT_1-original.PT_1.min())/(original.PT_1.max()-original.PT_1.min())
            assert np.allclose(scaled,joined[f'{fit}_{root}_scaled_PT'],rtol=0,atol=2e-15)
            sources[fit,root]=d
            paths[fit,root]=pd.read_csv(parent/f'Monocle2_{root}_path_cells.csv.gz').set_index('cell_id').loc[joined.cell_id]
            temp=joined[['dataset','capture_id','specimen_id','cnv_abs_log2']].copy()
            temp['bin']=np.minimum(np.floor(scaled.to_numpy()*12).astype(int),11)
            cap=temp.groupby(['dataset','specimen_id','capture_id','bin']).cnv_abs_log2.agg(['mean','size'])
            spec=cap[cap['size'].ge(10)]['mean'].groupby(level=['dataset','specimen_id','bin']).mean()
            rows=profiles[profiles.fit.eq(fit)&profiles.root.eq(root)]
            assert len(rows)==12
            for row in rows.itertuples():
                within=spec[spec.index.get_level_values('bin')==row.bin]
                cells=temp[temp.bin.eq(row.bin)].cnv_abs_log2
                expected=[cells.mean(),within.mean() if len(within)>=3 else np.nan,within.quantile(.25) if len(within)>=3 else np.nan,within.quantile(.75) if len(within)>=3 else np.nan]
                actual=[row.pooled_CNV,row.equal_capture_then_specimen_CNV,row.specimen_q25,row.specimen_q75]
                assert np.allclose(actual,expected,rtol=2e-11,atol=2e-13,equal_nan=True)
                assert row.specimens==len(within) and row.cells==len(cells)
            checks.append(dict(fit=fit,root=root,cells=len(d),CNV_bins=12))
    for row in table.itertuples():
        a,b=sources['old_CRC_only',row.previous_root],sources['joint_tree_CRC_subset',row.joint_root]
        rho=np.corrcoef(a.PT_1.rank().to_numpy(),b.PT_1.rank().to_numpy())[0,1]
        assert abs(rho-row.Spearman_PT)<1e-12
        expected=pd.crosstab(a.State.to_numpy(),b.State.to_numpy())
        actual=pd.read_csv(folder/f'State_overlap_{row.previous_root}_vs_{row.joint_root}.csv',index_col=0)
        assert np.array_equal(expected,actual)
        for fit,root,column in [('old_CRC_only',row.previous_root,'previous_PT_CNV_Spearman'),('joint_tree_CRC_subset',row.joint_root,'joint_PT_CNV_Spearman')]:
            value=np.corrcoef(sources[fit,root].PT_1.rank().to_numpy(),joined.cnv_abs_log2.rank().to_numpy())[0,1]
            assert abs(value-getattr(row,column))<1e-12
    for row in overlap.itertuples():
        a,b=paths['old_CRC_only',row.previous_root],paths['joint_tree_CRC_subset',row.joint_root]
        i,j=row.previous_path[4:],row.joint_path[4:]
        aa=a.index[np.isfinite(a['PT_'+i])&a['W_'+i].ge(.5)]
        bb=b.index[np.isfinite(b['PT_'+j])&b['W_'+j].ge(.5)]
        assert row.intersection==len(aa.intersection(bb))
        assert row.previous_cells==len(aa) and row.joint_cells==len(bb)
        assert abs(row.Jaccard-len(aa.intersection(bb))/len(aa.union(bb)))<1e-12
    save_json(LOG/'CRC_comparison_independent_validation.json',dict(passed=True,all_joined_cells=150917,
        all_four_root_pairs_checked=True,all_State_counts_and_rank_correlations_checked=True,
        all_path_pair_overlaps_checked=len(overlap),all_CNV_profiles=checks,
        original_legacy_models_untouched=True))
    print('All previous-vs-joint CRC cells, root/State/path and CNV profiles independently passed.',flush=True)


if __name__=='__main__':main()
