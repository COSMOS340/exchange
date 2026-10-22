"""Full-cell color-scale and library-size controls; no CNV value changes."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pipeline_common import ROOT,BLOCK,CELLS,TAB,LOG,save_json


def main():
    m=pd.read_csv(ROOT/'work'/BLOCK/'metadata.csv.gz',low_memory=False)
    cells=pd.read_csv(TAB/'Monocle2_primary_cells.csv.gz')
    assert np.array_equal(m.cell_id,cells.cell_id) and len(m)==CELLS
    xy=cells[['DDRTree_1','DDRTree_2']].to_numpy()
    centers=pd.read_csv(TAB/'Monocle2_centers.csv').set_index('center');centers.index=centers.index.astype(str)
    edges=pd.read_csv(TAB/'Monocle2_tree_edges.csv')
    cnv=m.cnv_abs_log2.to_numpy();lo,hi=np.quantile(cnv,[.01,.99])
    source=m[['cell_id','cnv_abs_log2','CNV_relative_band','nCount_input','Phase']].copy()
    source['log1p_nCount_input']=np.log1p(m.nCount_input)
    source.to_csv(TAB/'CNV_UMI_display_controls_all_cells.csv.gz',index=False)
    order=np.random.default_rng(20261022).permutation(CELLS)
    fig,axes=plt.subplots(2,2,figsize=(15,11),layout='constrained')
    for j,(ax,title,values) in enumerate(zip(axes.flat,
            ['Relative RNA CNV | original full range','Relative RNA CNV | 1st-99th percentile color limits',
             'Library size | log1p(nCount_input)','Inherited cell-cycle Phase'],
            [cnv,cnv,source.log1p_nCount_input.to_numpy(),m.Phase.to_numpy()])):
        if j<3:
            limits={'vmin':lo,'vmax':hi} if j==1 else {}
            sc=ax.scatter(*xy[order].T,c=values[order],s=.55,linewidths=0,cmap='viridis',rasterized=True,**limits)
            fig.colorbar(sc,ax=ax,shrink=.72,extend='both' if j==1 else 'neither',label='mean absolute log2' if j<2 else 'log1p input UMI')
        else:
            for level,color in [('G1','#909090'),('S','#E5A440'),('G2M','#607CB9')]:
                use=order[values[order]==level];ax.scatter(*xy[use].T,color=color,s=.55,linewidths=0,rasterized=True,label=f'{level}: {len(use):,}')
            ax.legend(frameon=False,fontsize=8)
        for a,b in edges[['from','to']].itertuples(index=False,name=None):
            line=centers.loc[[str(a),str(b)],['DDRTree_1','DDRTree_2']].to_numpy();ax.plot(*line.T,color='#333333',lw=.6,alpha=.7)
        span=np.ptp(xy,axis=0)
        ax.set(xlim=(xy[:,0].min()-.04*span[0],xy[:,0].max()+.04*span[0]),ylim=(xy[:,1].min()-.04*span[1],xy[:,1].max()+.04*span[1]),title=title,xlabel='DDRTree 1',ylabel='DDRTree 2')
        ax.set_aspect('equal',adjustable='box')
    fig.suptitle('One common full-cell tree | CNV, UMI and cycling controls\nOnly the top-right display colors are clipped; original CNV values and Low/Intermediate/High bands are unchanged.',fontsize=12)
    name='27_CNV_UMI_controls';out=ROOT/'results/01_图形/source'
    for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
    save_json(LOG/'CNV_UMI_control_figure_manifest.json',[dict(name=name,png='01_图形/source/'+name+'.png',pdf='01_图形/source/'+name+'.pdf',
        cells=CELLS,display_quantiles=[.01,.99],color_limits=[float(lo),float(hi)],
        below_display_limit=int((cnv<lo).sum()),above_display_limit=int((cnv>hi).sum()),
        original_CNV_values_and_bands_unchanged=True,sources=['02_表格/AllEpithelial/CNV_UMI_display_controls_all_cells.csv.gz'],visual_review='pending')])
    print('Saved full-cell CNV/UMI/Phase controls.',flush=True)


if __name__=='__main__':main()
