"""Display all preregistered K, root and cohort sensitivity results."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pipeline_common import ROOT,TAB,LOG,save_json


def main():
    folder=TAB/'dynamic_genes';sensitivity=folder/'sensitivity'
    out=ROOT/'results/01_图形/module_sensitivity';out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42})
    manifest=[];colors={'Module1':'#AC5B43','Module2':'#3C7F9C'}
    def save(fig,name,sources):
        for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=180,bbox_inches='tight')
        plt.close(fig);manifest.append(dict(name=name,png='01_图形/module_sensitivity/'+name+'.png',pdf='01_图形/module_sensitivity/'+name+'.pdf',sources=sources,visual_review='pending'))
    metrics=pd.read_csv(folder/'Monocle2_primary/dynamic_module_k_selection_metrics.csv')
    chosen=json.loads((folder/'Monocle2_primary/dynamic_module_selection.json').read_text())['chosen_k']
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for ax,field,title in zip(axes,['silhouette','median_bootstrap_ARI'],['Exact mean silhouette (every selected gene)','Median whole-linkage-cluster bootstrap ARI']):
        ax.plot(metrics.k,metrics[field],marker='o',color='#4D6E91')
        ax.axvline(chosen,color='#AC5B43',ls=':',label=f'Chosen K={chosen}')
        ax.set(xticks=range(2,9),xlabel='Candidate module count K',ylabel=field,title=title)
        ax.legend(frameon=False)
    axes[1].axhline(.75,color='#777777',ls='--',lw=.8)
    axes[1].set_ylim(0,1.03)
    fig.suptitle('Choose K before GO enrichment | maximum silhouette among stable candidates\n10 cohort-stratified linkage-cluster bootstrap runs; no cell sampling or GO-based tuning.',fontsize=11)
    save(fig,'32_module_K_stability',['02_表格/AllEpithelial/dynamic_genes/Monocle2_primary/dynamic_module_k_selection_metrics.csv','02_表格/AllEpithelial/dynamic_genes/Monocle2_primary/dynamic_module_bootstrap_stability.csv'])
    report=json.loads((LOG/'dynamic_sensitivity_summary.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for ax,root in zip(axes,['primary','alternative']):
        profile=pd.read_csv(sensitivity/f'{root}_fixed_primary_module_profiles.csv')
        for module,d in profile.groupby('module'):
            ax.plot((d.bin+.5)/30,d.mean_gene_Z,color=colors[module],label=module)
        ax.set(xlim=(0,1),xlabel='Within-root scaled pseudotime',ylabel='Mean member-gene Z of adjusted spline effect',title=root);ax.legend(frameon=False)
    fig.suptitle(f'Fixed primary gene memberships under both roots | same DDRTree\nSelected genes: {report["primary_selected"]:,} / {report["alternative_selected"]:,}; overlap {report["root_selection_overlap"]:,}; Jaccard {report["root_selection_Jaccard"]:.3f}. No root chosen by agreement.',fontsize=11)
    save(fig,'33_module_root_sensitivity',['02_表格/AllEpithelial/dynamic_genes/sensitivity/root_gene_selection_comparison.csv.gz',*[f'02_表格/AllEpithelial/dynamic_genes/sensitivity/{r}_fixed_primary_module_profiles.csv' for r in ['primary','alternative']]])
    summary=pd.read_csv(sensitivity/'all_cohort_leaveout_summary.csv')
    profile=pd.read_csv(sensitivity/'all_cohort_leaveout_fixed_module_profiles.csv')
    baseline=pd.read_csv(sensitivity/'primary_fixed_primary_module_profiles.csv')
    fig,axes=plt.subplots(2,4,figsize=(18,8),layout='constrained')
    for ax,row in zip(axes.flat,summary.itertuples()):
        if row.status!='assessed':
            ax.text(.5,.5,'Unassessable; retained failure record',ha='center',va='center',transform=ax.transAxes,wrap=True)
            ax.set_title('Leave out '+row.excluded_dataset);continue
        d=profile[profile.excluded_dataset.eq(row.excluded_dataset)]
        for module,v in d.groupby('module'):
            ax.plot((v.bin+.5)/30,v.mean_gene_Z,color=colors[module],label=module)
            b=baseline[baseline.module.eq(module)]
            ax.plot((b.bin+.5)/30,b.mean_gene_Z,color=colors[module],ls=':',alpha=.65)
        ax.set(xlim=(0,1),xlabel='Primary-root scaled pseudotime',ylabel='Mean member-gene Z',
               title=f'Leave out {row.excluded_dataset}\nGene Jaccard {row.selection_Jaccard:.2f}; partition ARI {row.primary_genes_partition_ARI:.2f}')
    fig.legend(handles=axes.flat[0].get_legend_handles_labels()[0],loc='outside lower center',ncol=2,frameon=False)
    fig.suptitle('All eight cohort leaveouts | refit within-specimen expression model, fixed common graph\nSolid: leaveout; dotted: full-pool reference. Fixed detection universe / primary memberships; no Harmony/tree/curve retraining.',fontsize=11)
    save(fig,'34_module_all_cohort_leaveout',['02_表格/AllEpithelial/dynamic_genes/sensitivity/all_cohort_leaveout_summary.csv','02_表格/AllEpithelial/dynamic_genes/sensitivity/all_cohort_leaveout_fixed_module_profiles.csv'])
    save_json(LOG/'module_root_cohort_figure_manifest.json',manifest)
    print('Saved K-selection, both roots and all eight cohort sensitivity figures.',flush=True)


if __name__=='__main__':main()
