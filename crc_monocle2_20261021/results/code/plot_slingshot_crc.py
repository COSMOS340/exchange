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
TAB=ROOT/'results/02_表格/CRC';OUT=ROOT/'results/01_图形/main';LOG=ROOT/'results/04_日志和检查'
m=pd.read_csv(ROOT/'work/CRC/metadata.csv.gz',low_memory=False)
assert len(m)==150917 and m.cell_id.is_unique
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
    ax.set(title=f'Slingshot {root} | same 150,917 CRC cells\n30D Harmony fit projected with embedCurves',xlabel='Inherited UMAP 1',ylabel='Inherited UMAP 2')
    ax.set_aspect('equal',adjustable='box')
    ax.spines[['top','right']].set_visible(False)
    name=f'12_Slingshot_{root}_UMAP'
    fig.savefig(OUT/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),dpi=180,bbox_inches='tight')
    plt.close(fig)
    assigned=data.filter(regex='^PT_').notna().any(axis=1)
    manifest.append(dict(name=name,png=f'01_图形/main/{name}.png',pdf=f'01_图形/main/{name}.pdf',
        question='同细胞Slingshot方法对照的路径及CNV分布如何？',
        finding=f'{assigned.sum():,}个细胞至少有一条有限拟时序，{(~assigned).sum():,}个未分配；全部150,917细胞保留显示。UMAP曲线为高维拟合后的投影，二维交叉不证明转换。',
        statistical_unit='all cells for geometry display; descriptive comparison, not independent validation',
        sources=[f'02_表格/CRC/Slingshot_{root}_cells.csv.gz',f'02_表格/CRC/Slingshot_{root}_curves.csv'],visual_review='pending'))
(LOG/'slingshot_figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
print('Saved Slingshot PNG/PDF pairs; visual review pending')
