"""Full-cell tissue plots and unit-equal source profiles on one saved tree."""
import json
import argparse
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from pipeline_common import ROOT,BLOCK,CELLS,TAB,LOG,save_json

SOURCE_COLORS={'Normal':'#379582','Adjacent':'#E5B04A','CRC':'#8060A8',
               'Normal_colorectal':'#379582','Primary_CRC':'#8060A8'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--only',nargs='*')
    args=parser.parse_args()
    m=pd.read_csv(ROOT/'work'/BLOCK/'metadata.csv.gz',low_memory=False)
    cells=pd.read_csv(TAB/'Monocle2_primary_cells.csv.gz')
    assert len(m)==CELLS and np.array_equal(m.cell_id,cells.cell_id)
    xy=cells[['DDRTree_1','DDRTree_2']].to_numpy()
    centers=pd.read_csv(TAB/'Monocle2_centers.csv').set_index('center')
    centers.index=centers.index.astype(str)
    edges=pd.read_csv(TAB/'Monocle2_tree_edges.csv')
    folder=TAB/'source_comparisons'
    out=ROOT/'results/01_图形/source';out.mkdir(exist_ok=True)
    order=np.random.default_rng(20261022).permutation(CELLS)
    span=np.ptp(xy,axis=0)
    limits=[(xy[:,j].min()-.04*span[j],xy[:,j].max()+.04*span[j]) for j in range(2)]
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    manifest=[]
    def tree(ax):
        for a,b in edges[['from','to']].itertuples(index=False,name=None):
            line=centers.loc[[str(a),str(b)],['DDRTree_1','DDRTree_2']].to_numpy()
            ax.plot(*line.T,color='#292929',lw=.65,alpha=.65,zorder=3)
    def points(ax,field,title,selected=None):
        valid=np.ones(CELLS,dtype=bool) if selected is None else np.asarray(selected)
        if selected is not None:
            ax.scatter(*xy.T,c='#EAEAEA',s=.3,linewidths=0,rasterized=True,marker='o')
        use=order[valid[order]]
        ax.scatter(*xy[use].T,c=m[field].map(SOURCE_COLORS).to_numpy()[use],s=.7,
                   linewidths=0,alpha=.75,rasterized=True,marker='o')
        tree(ax);ax.set(xlim=limits[0],ylim=limits[1],title=title,xlabel='DDRTree 1',ylabel='DDRTree 2')
        ax.set_aspect('equal',adjustable='box')
    def legend(fig,levels):
        handles=[Line2D([],[],marker='o',ls='',color=SOURCE_COLORS[x],
                       label=x,markersize=5) for x in levels]
        fig.legend(handles=handles,loc='outside lower center',ncol=len(levels),frameon=False)
    def save(fig,name,question,sources,unit):
        fig.canvas.draw();fig.set_layout_engine('none')
        if not args.only or name in args.only:
            fig.savefig(out/(name+'.png'),dpi=180,bbox_inches='tight')
            fig.savefig(out/(name+'.pdf'),dpi=180,bbox_inches='tight')
        else:
            assert (out/(name+'.png')).exists() and (out/(name+'.pdf')).exists()
        plt.close(fig)
        manifest.append(dict(name=name,png='01_图形/source/'+name+'.png',pdf='01_图形/source/'+name+'.pdf',
            question=question,sources=sources,statistical_unit=unit,fitted_cells=CELLS,
            finding='Pending actual figure inspection',visual_review='pending'))
    for number,field,levels in [(20,'tissue_group',['Normal_colorectal','Primary_CRC']),
                                (21,'display_3class',['Normal','Adjacent','CRC'])]:
        fig,axes=plt.subplots(1,len(levels)+1,figsize=(5*(len(levels)+1),4.7),layout='constrained')
        points(axes[0],field,'All 223,736 epithelial cells')
        for ax,level in zip(axes[1:],levels):
            keep=m[field].eq(level)
            points(ax,field,level+' | n='+format(int(keep.sum()),','),keep)
        legend(fig,levels)
        subtitle=('Normal_colorectal includes unspecified-distance normal and explicitly adjacent tissue; no healthy-donor category' if number==20 else 'Normal = unspecified-distance normal; Adjacent = explicitly adjacent; no healthy-donor category')
        fig.suptitle('Frozen tissue sources | one common expression tree\n'+subtitle,fontsize=11)
        save(fig,str(number)+'_tissue_source_facets','组织来源在完全相同的共同树上如何分布？',
             ['02_表格/AllEpithelial/source_mapping_all_cells.csv.gz','02_表格/AllEpithelial/Monocle2_primary_cells.csv.gz'],'All cells; display only, no tissue-based refit')
    cohorts=sorted(m.dataset.unique())
    fig,axes=plt.subplots(2,4,figsize=(19,9),layout='constrained')
    for ax,cohort in zip(axes.flat,cohorts):
        keep=m.dataset.eq(cohort);points(ax,'display_3class',cohort+' | n='+format(int(keep.sum()),','),keep)
    legend(fig,['Normal','Adjacent','CRC'])
    fig.suptitle('Cohort composition on the same full tree | Adjacent occurs only in GSE245552',fontsize=12)
    save(fig,'22_source_by_cohort','各来源区段是否与队列构成重合？',
         ['02_表格/AllEpithelial/source_counts_by_dataset.csv'],'Full-cell cohort facets, same coordinates and tree')
    coverage=pd.read_csv(folder/'shared_identity_source_coverage.csv')
    shared=coverage.loc[coverage.shared_display_eligible,'broad_identity'].tolist()
    rows=max(1,int(np.ceil(len(shared)/3)))
    fig,axes=plt.subplots(rows,3,figsize=(18,4.8*rows),layout='constrained',squeeze=False)
    for ax,identity in zip(axes.flat,shared):
        keep=m.broad_identity.eq(identity)
        points(ax,'display_3class',textwrap.fill(identity,32)+' | n='+format(int(keep.sum()),','),keep)
    for ax in list(axes.flat)[len(shared):]:ax.axis('off')
    if not shared:
        axes.flat[0].text(.5,.5,'No identity meets the stated shared-source coverage thresholds',
                          ha='center',va='center',transform=axes.flat[0].transAxes,wrap=True)
    legend(fig,['Normal','Adjacent','CRC'])
    fig.suptitle('Shared broad identities | >=30 cells and >=3 captures in each two-class source',fontsize=12)
    save(fig,'23_shared_identity_sources','共享宽上皮身份内，组织来源分布是否仍然不同？',
         ['02_表格/AllEpithelial/source_comparisons/shared_identity_source_coverage.csv'],'Identity subsets on full-cell background; no refit or new eligibility filtering')
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    for row,root in enumerate(['primary','alternative']):
        table=pd.read_csv(folder/(root+'_source_PT_distributions.csv'))
        for source,d in table.groupby('source',sort=False):
            for col,feature in enumerate(['pooled_mass','equal_capture_then_specimen_mass']):
                axes[row,col].plot((d.bin+.5)/30,d[feature],color=SOURCE_COLORS[source],label=source,lw=1.5)
        for col in range(2):axes[row,col].set(xlabel='Within-root scaled pseudotime',ylabel='Within-source fraction per bin',title=root+' | '+('pooled cells' if col==0 else 'equal capture -> specimen'))
    legend(fig,['Normal','Adjacent','CRC'])
    save(fig,'24_source_time_distributions','来源的拟时序分布在细胞合并与标本等权后是否相同？',
         ['02_表格/AllEpithelial/source_comparisons/'+r+'_source_PT_distributions.csv' for r in ['primary','alternative']],
         'True bin probability masses; equal captures within specimens, then equal specimens within source')
    identities=['All_identities',*shared]
    fig,axes=plt.subplots(len(identities),2,figsize=(14,3.1*len(identities)),layout='constrained',squeeze=False)
    for col,root in enumerate(['primary','alternative']):
        profiles=pd.read_csv(folder/(root+'_source_CNV_profiles.csv.gz'))
        for row,identity in enumerate(identities):
            ax=axes[row,col]
            for source,d in profiles[profiles.identity.eq(identity)].groupby('source',sort=False):
                x=(d.bin.to_numpy()+.5)/30
                ax.plot(x,d['mean'],color=SOURCE_COLORS[source],lw=1.3,marker='o',markersize=2.3,label=source)
                ax.fill_between(x,d.q25,d.q75,color=SOURCE_COLORS[source],alpha=.15)
            ax.set(title=root+' | '+textwrap.fill(identity,36),xlabel='Scaled pseudotime',ylabel='Relative RNA CNV\nmean absolute log2',xlim=(0,1))
    for row in range(len(identities)):
        bounds=[ax.get_ylim() for ax in axes[row]]
        common=(min(x[0] for x in bounds),max(x[1] for x in bounds))
        for ax in axes[row]:ax.set_ylim(common)
    legend(fig,['Normal','Adjacent','CRC'])
    fig.suptitle('Equal capture -> specimen CNV profiles | mean / between-specimen IQR\nBins require >=3 specimens; Adjacent is confounded with GSE245552. IQR is not a confidence interval.',fontsize=11)
    save(fig,'26_source_CNV_within_shared_identities','共同轨迹的CNV差异在共享身份和来源内是否仍出现？',
         ['02_表格/AllEpithelial/source_comparisons/'+r+'_source_CNV_profiles.csv.gz' for r in ['primary','alternative']],
         'Capture -> specimen equal mean/IQR, at least3 specimens per bin; descriptive, no clinical time or causal claim')
    for root in ['primary','alternative']:
        table=pd.read_csv(folder/(root+'_leave_one_cohort_source_profiles.csv'))
        fig,axes=plt.subplots(2,4,figsize=(19,9),layout='constrained')
        for ax,excluded in zip(axes.flat,cohorts):
            for source,d in table[table.excluded_dataset.eq(excluded)].groupby('source',sort=False):
                ax.plot((d.bin+.5)/30,d.equal_specimen_fraction,color=SOURCE_COLORS[source],label=source)
            ax.set(title='Leave out '+excluded,xlabel='Scaled pseudotime',ylabel='Equal-specimen source fraction',ylim=(-.02,1.02))
        legend(fig,['Normal','Adjacent','CRC'])
        fig.suptitle(root+' | fixed-tree descriptive leaveout; no Harmony/curve/tree refit\nDropping GSE245552 removes all Adjacent observations',fontsize=12)
        save(fig,'25_'+root+'_source_cohort_leaveout','逐一去除队列后来源组成描述是否改变？',
             ['02_表格/AllEpithelial/source_comparisons/'+root+'_leave_one_cohort_source_profiles.csv'],
             'Capture -> specimen equal fractions; insufficient bins remain gaps, no inferential test')
    save_json(LOG/'source_figure_manifest.json',manifest)
    print('Saved',len(manifest),'source PNG/PDF pairs; actual visual review still pending.',flush=True)


if __name__=='__main__':main()
