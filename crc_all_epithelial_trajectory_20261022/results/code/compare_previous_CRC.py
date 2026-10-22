"""Exact-cell comparison with frozen CRC-only results, without refitting them."""
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from pipeline_common import ROOT,BLOCK,TAB,LOG,digest,save_json


def scale(x):return (x-x.min())/(x.max()-x.min())


def main():
    old=ROOT.parent.parent/'exchange/crc_monocle2_20261021/results/02_表格/CRC'
    out=TAB/'CRC_only_comparison';out.mkdir(exist_ok=False)
    sources=out/'frozen_previous_outputs';sources.mkdir()
    m=pd.read_csv(ROOT/'work'/BLOCK/'metadata.csv.gz',low_memory=False)
    crc=m[m.tissue_group.eq('Primary_CRC')].reset_index(drop=True)
    assert len(crc)==150917 and crc.cell_id.is_unique
    joined=crc[['cell_id','dataset','capture_id','specimen_id','subject_key','reported_patient_link_supported','broad_identity','cnv_abs_log2','CNV_relative_band']].copy()
    records=[];tables={};paths={};stats=[];profiles=[];path_rows=[]
    for fit in ['old_CRC_only','joint_tree_CRC_subset']:
        for root in ['primary','alternative']:
            source=(old if fit=='old_CRC_only' else TAB)/f'Monocle2_{root}_cells.csv.gz'
            d=pd.read_csv(source).set_index('cell_id')
            assert set(crc.cell_id)<=set(d.index)
            d=d.loc[crc.cell_id].reset_index();assert np.array_equal(d.cell_id,crc.cell_id)
            if fit=='old_CRC_only':
                assert len(d)==len(pd.read_csv(source))
                shutil.copy2(source,sources/source.name)
                records.append(dict(path='02_表格/AllEpithelial/CRC_only_comparison/frozen_previous_outputs/'+source.name,bytes=source.stat().st_size,sha256=digest(source)))
            tables[fit,root]=d
            joined[f'{fit}_{root}_State']=d.State.to_numpy()
            joined[f'{fit}_{root}_PT']=d.PT_1.to_numpy()
            # For the joint fit, scale by the ORIGINAL FULL-POOL distance range.
            # CRC subset need not span the endpoint extremes.
            full=pd.read_csv(source,usecols=['PT_1']).PT_1
            pt=(d.PT_1-full.min())/(full.max()-full.min())
            joined[f'{fit}_{root}_scaled_PT']=pt.to_numpy()
            temp=crc[['dataset','capture_id','specimen_id','cnv_abs_log2']].copy()
            temp['bin']=np.minimum((pt.to_numpy()*12).astype(int),11)
            capture=temp.groupby(['dataset','specimen_id','capture_id','bin']).agg(cells=('cnv_abs_log2','size'),value=('cnv_abs_log2','mean')).reset_index()
            capture['eligible']=capture.cells.ge(10)
            capture.to_csv(out/f'{fit}_{root}_CNV_capture_bins.csv',index=False)
            spec=capture[capture.eligible].groupby(['dataset','specimen_id','bin']).agg(value=('value','mean'),captures=('capture_id','nunique'),cells=('cells','sum')).reset_index()
            spec.to_csv(out/f'{fit}_{root}_CNV_specimen_bins.csv',index=False)
            for b in range(12):
                s=spec[spec.bin.eq(b)];p=temp[temp.bin.eq(b)]
                profiles.append(dict(fit=fit,root=root,bin=b,cells=len(p),specimens=s.specimen_id.nunique(),captures=int(s.captures.sum()),
                    pooled_CNV=float(p.cnv_abs_log2.mean()) if len(p) else np.nan,
                    equal_capture_then_specimen_CNV=float(s.value.mean()) if len(s)>=3 else np.nan,
                    specimen_q25=float(s.value.quantile(.25)) if len(s)>=3 else np.nan,
                    specimen_q75=float(s.value.quantile(.75)) if len(s)>=3 else np.nan))
            psource=(old if fit=='old_CRC_only' else TAB)/f'Monocle2_{root}_path_cells.csv.gz'
            path=pd.read_csv(psource).set_index('cell_id').loc[crc.cell_id].reset_index()
            paths[fit,root]=path
            if fit=='old_CRC_only':
                shutil.copy2(psource,sources/psource.name)
                records.append(dict(path='02_表格/AllEpithelial/CRC_only_comparison/frozen_previous_outputs/'+psource.name,bytes=psource.stat().st_size,sha256=digest(psource)))
    for oldroot in ['primary','alternative']:
        for newroot in ['primary','alternative']:
            a,b=tables['old_CRC_only',oldroot],tables['joint_tree_CRC_subset',newroot]
            stats.append(dict(previous_root=oldroot,joint_root=newroot,cells=len(a),Spearman_PT=float(spearmanr(a.PT_1,b.PT_1).statistic),
                previous_PT_CNV_Spearman=float(spearmanr(a.PT_1,crc.cnv_abs_log2).statistic),
                joint_PT_CNV_Spearman=float(spearmanr(b.PT_1,crc.cnv_abs_log2).statistic)))
            pd.crosstab(a.State,b.State).to_csv(out/f'State_overlap_{oldroot}_vs_{newroot}.csv')
            pa,pb=paths['old_CRC_only',oldroot],paths['joint_tree_CRC_subset',newroot]
            for ca in [x for x in pa if x.startswith('PT_')]:
                ia=ca.split('_')[1];ma=np.isfinite(pa[ca])&pa['W_'+ia].ge(.5)
                for cb in [x for x in pb if x.startswith('PT_')]:
                    ib=cb.split('_')[1];mb=np.isfinite(pb[cb])&pb['W_'+ib].ge(.5)
                    path_rows.append(dict(previous_root=oldroot,joint_root=newroot,previous_path='Path'+ia,joint_path='Path'+ib,
                        previous_cells=int(ma.sum()),joint_cells=int(mb.sum()),intersection=int((ma&mb).sum()),Jaccard=float((ma&mb).sum()/(ma|mb).sum())))
    joined.to_csv(out/'exact_150917_cell_comparison.csv.gz',index=False)
    pd.DataFrame(stats).to_csv(out/'all_four_root_comparisons.csv',index=False)
    pd.DataFrame(path_rows).to_csv(out/'all_path_pair_overlaps.csv',index=False)
    profile=pd.DataFrame(profiles);profile.to_csv(out/'all_CNV_binned_profiles.csv',index=False)
    save_json(LOG/'CRC_only_comparison_validation.json',dict(passed=True,cells=150917,exact_one_to_one_cell_IDs=True,
        previous_models_refitted=False,previous_release='https://github.com/COSMOS340/exchange/releases/tag/crc-monocle2-pc-results-v3-20261022',
        previous_result_code_commit='4d8a13008e57dc5a30a6fb5b5fc1a48410151cce',frozen_previous_table_manifest=records,
        all_four_root_pairs=stats,path_memberships_overlap_on_shared_trunks=True,fit_specific_State_labels=True))
    plot_comparisons(joined,stats,profile,out)
    print('All exact-cell old/new root, State, path and CNV comparisons saved.',flush=True)


def plot_comparisons(joined,stats,profile,out):
    figout=ROOT/'results/01_图形/CRC_only_comparison';figout.mkdir(exist_ok=True)
    plt.rcParams.update({'pdf.fonttype':42,'font.size':9})
    manifest=[]
    def save(fig,name):
        for ext in ['png','pdf']:fig.savefig(figout/(name+'.'+ext),dpi=180,bbox_inches='tight')
        plt.close(fig);manifest.append(dict(name=name,png='01_图形/CRC_only_comparison/'+name+'.png',pdf='01_图形/CRC_only_comparison/'+name+'.pdf',
            sources=['02_表格/AllEpithelial/CRC_only_comparison/exact_150917_cell_comparison.csv.gz','02_表格/AllEpithelial/CRC_only_comparison/all_CNV_binned_profiles.csv'],visual_review='pending'))
    fig,axes=plt.subplots(2,2,figsize=(11,10),layout='constrained')
    for ax,stat in zip(axes.flat,stats):
        artist=ax.hexbin(joined[f'old_CRC_only_{stat["previous_root"]}_scaled_PT'],joined[f'joint_tree_CRC_subset_{stat["joint_root"]}_scaled_PT'],gridsize=80,mincnt=1,bins='log',cmap='viridis')
        fig.colorbar(artist,ax=ax,label='Cell count (log scale)',shrink=.75)
        ax.set(xlabel='Frozen CRC-only scaled PT',ylabel='Joint-tree CRC-subset scaled PT',title=f'old {stat["previous_root"]} / joint {stat["joint_root"]} | rho={stat["Spearman_PT"]:.3f}',xlim=(0,1),ylim=(0,1))
    fig.suptitle('Exact same 150,917 CRC cell IDs | all four root pairs\nHexagons aggregate every cell; density shown on a log scale. No root-pair selection or p-values.')
    save(fig,'40_CRC_only_all_root_pairs')
    fig,ax=plt.subplots(figsize=(10,8),layout='constrained')
    ct=pd.read_csv(out/'State_overlap_primary_vs_primary.csv',index_col=0)
    for stat in stats:
        other=pd.read_csv(out/f'State_overlap_{stat["previous_root"]}_vs_{stat["joint_root"]}.csv',index_col=0)
        assert ct.equals(other), 'Native State assignments differ across rerooting; cannot collapse panels'
    normalized=ct.div(ct.sum(axis=1),axis=0)
    im=ax.imshow(normalized,cmap='Blues',vmin=0,vmax=1,aspect='auto')
    ax.set(xticks=range(len(ct.columns)),xticklabels=ct.columns,yticks=range(len(ct.index)),yticklabels=ct.index,
           xlabel='Joint-tree State (same CRC subset)',ylabel='Frozen CRC-only State')
    for i in range(len(ct)):
        for j in range(len(ct.columns)):
            ax.text(j,i,f'{ct.iloc[i,j]:,}',ha='center',va='center',fontsize=10,color='white' if normalized.iloc[i,j]>.6 else 'black')
    fig.colorbar(im,ax=ax,label='Within-old-State fraction',shrink=.8)
    fig.suptitle('State IDs are fit-specific | exact same 150,917 CRC cells\nNative State assignments are invariant under rerooting within each fit; all four root-pair tables are identical.',fontsize=11)
    save(fig,'41_CRC_only_State_overlap')
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for ax,root in zip(axes,['primary','alternative']):
        for fit,color in [('old_CRC_only','#CC8D31'),('joint_tree_CRC_subset','#457CB7')]:
            d=profile[profile.root.eq(root)&profile.fit.eq(fit)];x=(d.bin.to_numpy()+.5)/12
            ax.plot(x,d.equal_capture_then_specimen_CNV,label=fit,color=color)
            ax.fill_between(x,d.specimen_q25,d.specimen_q75,color=color,alpha=.12)
            ax.plot(x,d.pooled_CNV,color=color,ls=':',alpha=.7)
        ax.set(xlabel='Within-fit full-pool scaled PT',ylabel='Relative RNA CNV (mean absolute log2)',title=root)
    axes[1].legend(fontsize=8);fig.suptitle('Same frozen CNV values; different fitted geometry/order\nSolid: equal capture -> specimen mean / IQR (>=3 specimens). Dotted: pooled cells. IQR is not a confidence interval.',fontsize=11)
    save(fig,'42_CRC_only_CNV_profiles')
    save_json(LOG/'CRC_comparison_figure_manifest.json',manifest)


if __name__=='__main__':
    import sys
    if '--plots-only' in sys.argv:
        out=TAB/'CRC_only_comparison'
        plot_comparisons(pd.read_csv(out/'exact_150917_cell_comparison.csv.gz'),pd.read_csv(out/'all_four_root_comparisons.csv').to_dict('records'),pd.read_csv(out/'all_CNV_binned_profiles.csv'),out)
    else:main()
