"""Plot saved 30D Slingshot curves after their native UMAP projection."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
TAB=ROOT/'results/02_表格/AllEpithelial';OUT=ROOT/'results/01_图形/main';LOG=ROOT/'results/04_日志和检查'
m=pd.read_csv(ROOT/'work/AllEpithelial/metadata.csv.gz',low_memory=False)
assert len(m)==223736 and m.cell_id.is_unique
xy=m[['state_UMAP_1','state_UMAP_2']].to_numpy()
assert np.isfinite(xy).all()
palette={'Low':'#3182BD','Intermediate':'#BDBDBD','High':'#D94B4B'}
order=np.random.default_rng(20261020).permutation(len(m))
manifest=[]
for root in ['primary','alternative']:
    data=pd.read_csv(TAB/f'Slingshot_{root}_cells.csv.gz')
    assert data.cell_id.is_unique and set(data.cell_id)==set(m.cell_id)
    data=data.set_index('cell_id').loc[m.cell_id]
    curves=pd.read_csv(TAB/f'Slingshot_{root}_curves.csv')
    fig,ax=plt.subplots(figsize=(11,8),layout='constrained')
    ax.scatter(*xy[order].T,c=m.CNV_relative_band.map(palette).to_numpy()[order],
               s=.7,marker='o',linewidths=0,rasterized=True)
    for number,curve in curves.groupby('path'):
        ax.plot(curve.UMAP_1,curve.UMAP_2,lw=1.1,label=f'Path {number}')
    handles,labels=ax.get_legend_handles_labels()
    handles += [Line2D([],[],marker='o',linestyle='',color=c,label=f'{b}: n={m.CNV_relative_band.eq(b).sum():,}') for b,c in palette.items()]
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.01,1),frameon=False)
    ax.set(title=f'Slingshot {root} | same 223,736 epithelial cells\n30D Harmony fit projected with embedCurves',xlabel='Inherited UMAP 1',ylabel='Inherited UMAP 2')
    ax.set_aspect('equal',adjustable='box')
    ax.spines[['top','right']].set_visible(False)
    name=f'12_Slingshot_{root}_UMAP'
    fig.savefig(OUT/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),dpi=180,bbox_inches='tight')
    plt.close(fig)
    assigned=data.filter(regex='^PT_').notna().any(axis=1)
    manifest.append(dict(name=name,png=f'01_图形/main/{name}.png',pdf=f'01_图形/main/{name}.pdf',
        question='同细胞Slingshot方法对照的路径及CNV分布如何？',
        finding=f'{assigned.sum():,}个细胞至少有一条有限拟时序，{(~assigned).sum():,}个未分配；全部223,736细胞保留显示。UMAP曲线为高维拟合后的投影，二维交叉不证明转换。',
        statistical_unit='all cells for geometry display; descriptive comparison, not independent validation',
        sources=[f'02_表格/AllEpithelial/Slingshot_{root}_cells.csv.gz',f'02_表格/AllEpithelial/Slingshot_{root}_curves.csv'],visual_review='pending'))
    fig,axes=plt.subplots(1,4,figsize=(20,5.2),layout='constrained')
    sources={'Normal':'#379582','Adjacent':'#E5B04A','CRC':'#8060A8'}
    span=np.ptp(xy,axis=0)
    for ax,level in zip(axes,['All','Normal','Adjacent','CRC']):
        selected=np.ones(len(m),bool) if level=='All' else m.display_3class.eq(level).to_numpy()
        ax.scatter(*xy.T,c='#E8E8E8',s=.3,linewidths=0,rasterized=True)
        use=order[selected[order]]
        ax.scatter(*xy[use].T,c=m.display_3class.map(sources).to_numpy()[use],s=.65,linewidths=0,rasterized=True)
        for _,curve in curves.groupby('path'):ax.plot(curve.UMAP_1,curve.UMAP_2,color='#222222',lw=.6,alpha=.55)
        ax.set(title=f'{level} | n={int(selected.sum()):,}',xlabel='Inherited UMAP 1',ylabel='Inherited UMAP 2',
               xlim=(xy[:,0].min()-.04*span[0],xy[:,0].max()+.04*span[0]),ylim=(xy[:,1].min()-.04*span[1],xy[:,1].max()+.04*span[1]))
        ax.set_aspect('equal',adjustable='box')
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color=c,label=x) for x,c in sources.items()],loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle(f'Slingshot {root} | same full-cell fit / same native UMAP projection in every facet\nNormal: unspecified-distance normal; Adjacent: explicitly adjacent (only GSE245552); no healthy-donor category.',fontsize=11)
    name=f'13_Slingshot_{root}_tissue_sources'
    for ext in ['png','pdf']:fig.savefig(OUT/(name+'.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
    manifest.append(dict(name=name,png=f'01_图形/main/{name}.png',pdf=f'01_图形/main/{name}.pdf',
        sources=['02_表格/AllEpithelial/source_mapping_all_cells.csv.gz',f'02_表格/AllEpithelial/Slingshot_{root}_curves.csv'],
        statistical_unit='All 223736 cells displayed; sources are annotations, not separate fits',visual_review='pending'))
(LOG/'slingshot_figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
print('Saved Slingshot PNG/PDF pairs; visual review pending')
