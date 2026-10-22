"""Complete selected-gene heatmaps, frozen annotations and GO source links."""
import textwrap
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from pipeline_common import ROOT,TAB,OBJ,LOG,save_json


def main():
    folder=TAB/'dynamic_genes/Monocle2_primary'
    assignments=pd.read_csv(folder/'dynamic_gene_module_assignment.csv')
    reps=pd.read_csv(folder/'dynamic_module_representative_genes.csv')
    go=pd.read_csv(folder/'dynamic_module_GO_BP_display.csv').fillna('')
    annotation=pd.read_csv(TAB/'source_comparisons/primary_dynamic_heatmap_annotations.csv.gz')
    coverage=pd.read_csv(folder/'dynamic_bin_coverage.csv')
    with h5py.File(OBJ/'dynamic_model_primary.h5') as h:
        genes=h['gene'].asstr()[:];z=h['gene_Z'][:].T
    assert np.array_equal(genes,assignments.gene)
    selected=assignments[assignments.selected].copy()
    selected['peak_bin']=[int(np.nanargmax(z[i])) for i in selected.index]
    selected=selected.sort_values(['module','peak_bin','gene'])
    selected[['gene','module','peak_bin']].to_csv(folder/'dynamic_heatmap_gene_order.csv',index=False)
    fields=['tissue_group','display_3class','CNV_relative_band','broad_identity','dataset']
    colors={};pal=plt.get_cmap('tab20')
    fixed={'Normal_colorectal':'#379582','Primary_CRC':'#8060A8','Normal':'#379582',
           'Adjacent':'#E5B04A','CRC':'#8060A8','Low':'#278a4f','Intermediate':'#E4C24A','High':'#A957BB'}
    for field in fields:
        levels=annotation.loc[annotation.field.eq(field),'level'].drop_duplicates()
        colors[field]={x:fixed.get(x,pal(i)) for i,x in enumerate(levels)}
    out=ROOT/'results/01_图形/dynamic_genes';out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42})
    manifest=[]
    for name,rows in [('30_dynamic_gene_heatmap',selected)]+[
            ('31_'+module+'_zoom',selected[selected.module.eq(module)]) for module in sorted(selected.module.unique())]:
        fig=plt.figure(figsize=(18,14),layout='constrained')
        grid=fig.add_gridspec(7,2,width_ratios=[11,6],height_ratios=[.42]*5+[8,.45])
        for j,field in enumerate(fields):
            ax=fig.add_subplot(grid[j,0]);bottom=np.zeros(30)
            d=annotation[annotation.field.eq(field)]
            for level,color in colors[field].items():
                values=d[d.level.eq(level)].set_index('bin').equal_capture_then_specimen_fraction.reindex(range(30)).to_numpy()
                ax.bar(np.arange(30)+.5,values,bottom=bottom,width=1,color=color,linewidth=0)
                bottom+=values
            ax.set(xlim=(0,30),ylim=(0,1),ylabel=['2-class','3-class','CNV','Identity','Cohort'][j])
            ax.yaxis.label.set(rotation=0,ha='right',va='center');ax.set_xticks([]);ax.set_yticks([])
        for j,field in enumerate(fields):
            legend_ax=fig.add_subplot(grid[j,1]);legend_ax.axis('off')
            handles=[Patch(facecolor=c,label=textwrap.fill(x,22)) for x,c in colors[field].items()]
            legend_ax.legend(handles=handles,loc='center left',ncol=3 if len(handles)>3 else len(handles),
                             frameon=False,fontsize=6.5,handlelength=1,columnspacing=.8)
        ax=fig.add_subplot(grid[5,0]);matrix=z[rows.index]
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#D6D6D6')
        im=ax.imshow(matrix,aspect='auto',interpolation='nearest',extent=(0,30,len(rows),0),cmap=cmap,norm=Normalize(-2.5,2.5))
        boundaries=[];offset=0
        for module,d in rows.groupby('module',sort=False):
            if offset:ax.axhline(offset,color='white',lw=1.4)
            boundaries.append((offset,offset+len(d),module));offset+=len(d)
        ax.set(xlim=(0,30),xlabel='Primary-root scaled pseudotime (30 fixed bins)',ylabel='All selected genes; row Z of adjusted spline effect')
        ax.set_xticks([0,6,12,18,24,30],['0','.2','.4','.6','.8','1'])
        ax.set_yticks([(a+b)/2 for a,b,_ in boundaries],[f'{m} | {b-a} genes' for a,b,m in boundaries])
        sidebar=fig.add_subplot(grid[5,1]);sidebar.axis('off')
        for a,b,module in boundaries:
            midpoint=1-(a+b)/2/len(rows)
            representative=', '.join(reps.loc[reps.module.eq(module)].sort_values('rank').gene)
            lines=[f'{module}: {b-a} genes',f'Centroid-nearest genes: {representative}','GO Biological Process:']
            for item in go[go.module.eq(module)].itertuples():
                line=(f'{item.GO_ID} {item.term} | q={item.BH_q_value:.2g}' if item.GO_ID else item.term)
                lines.extend(textwrap.wrap(line,53))
            sidebar.text(.03,midpoint,'\n'.join(lines),va='center',fontsize=9,linespacing=1.5)
        bar=fig.add_subplot(grid[6,0]);fig.colorbar(im,cax=bar,orientation='horizontal',extend='both',label='Within-gene Z; display clipped at +/-2.5 (original values retained)')
        note=fig.add_subplot(grid[6,1]);note.axis('off');note.text(0,.5,'GO: measured detection-eligible background, BH across module-term tests.\nAssociation conditional on the inferred tree; no real-time or causal interpretation.',va='center',fontsize=8)
        fig.suptitle(f'All epithelial cells | {len(rows):,} dynamic genes\nSpecimen-intercept-adjusted log1p(CP10K); q <= .05, amplitude >= .25\nTop strips: equal capture -> specimen fractions; gray column = insufficient coverage',fontsize=12)
        for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=180,bbox_inches='tight')
        plt.close(fig)
        manifest.append(dict(name=name,png='01_图形/dynamic_genes/'+name+'.png',pdf='01_图形/dynamic_genes/'+name+'.pdf',
            genes=len(rows),finite_values=int(np.isfinite(matrix).sum()),color_clipped_values=int((np.abs(matrix)>2.5).sum()),
            unsupported_bins=coverage.loc[~coverage.supported,'bin'].tolist(),
            sources=['02_表格/AllEpithelial/dynamic_genes/Monocle2_primary/dynamic_gene_smoothed_expression.csv.gz',
                     '02_表格/AllEpithelial/dynamic_genes/Monocle2_primary/dynamic_heatmap_gene_order.csv',
                     '02_表格/AllEpithelial/source_comparisons/primary_dynamic_heatmap_annotations.csv.gz',
                     '02_表格/AllEpithelial/dynamic_genes/Monocle2_primary/dynamic_module_GO_BP_full.csv'],
            visual_review='pending'))
    pd.DataFrame([dict(field=f,level=x,color=matplotlib.colors.to_hex(c)) for f,d in colors.items() for x,c in d.items()]).to_csv(folder/'dynamic_heatmap_annotation_colors.csv',index=False)
    save_json(LOG/'dynamic_heatmap_figure_manifest.json',manifest)
    print('Saved full gene heatmap and both complete module zooms.',flush=True)


if __name__=='__main__':main()
