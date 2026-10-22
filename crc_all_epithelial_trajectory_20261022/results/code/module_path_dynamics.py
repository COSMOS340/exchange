"""All-cell fixed module activities and representative genes on every path."""
import json
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pipeline_common import ROOT,BLOCK,CELLS,TAB,OBJ,LOG,save_json


def calculate_cells():
    dest=TAB/'dynamic_genes/module_and_representative_values_all_cells.csv.gz'
    if dest.exists():raise RuntimeError('All-cell module values already exist; validate/version them')
    m=pd.read_csv(ROOT/'work'/BLOCK/'metadata.csv.gz',low_memory=False)
    genes=pd.read_csv(ROOT/'work'/BLOCK/'genes.csv').gene.to_numpy()
    standard=pd.read_csv(TAB/'dynamic_genes/all_gene_expression_standardization.csv')
    assert np.array_equal(genes,standard.gene)
    assignment=pd.read_csv(TAB/'dynamic_genes/Monocle2_primary/dynamic_gene_module_assignment.csv')
    reps=pd.read_csv(TAB/'dynamic_genes/Monocle2_primary/dynamic_module_representative_genes.csv')
    modules=sorted(assignment.loc[assignment.selected,'module'].unique())
    feature_names=modules+reps.sort_values(['module','rank']).gene.tolist()
    weights=np.zeros((len(genes),len(feature_names)));offset=np.zeros(len(feature_names))
    for j,module in enumerate(modules):
        members=assignment.module.eq(module).to_numpy();assert members.sum()>0
        sd=standard.population_sd.to_numpy()[members];assert (sd>1e-12).all()
        weights[members,j]=1/sd/members.sum()
        offset[j]=-np.mean(standard.mean_log1p_CP10K.to_numpy()[members]/sd)
    for j,gene in enumerate(feature_names[len(modules):],len(modules)):
        weights[np.flatnonzero(genes==gene),j]=1
    values=np.empty((CELLS,len(feature_names)))
    with h5py.File(ROOT/'work'/BLOCK/'counts.h5') as h:
        ptr=h['indptr'][:]
        for lo in range(0,CELLS,3001):
            hi=min(lo+3001,CELLS);a,b=int(ptr[lo]),int(ptr[hi])
            x=sp.csr_matrix((h['data'][a:b].astype(float),h['indices'][a:b],ptr[lo:hi+1]-a),shape=(hi-lo,len(genes)))
            x=x.multiply((1e4/m.nCount_input.iloc[lo:hi].to_numpy())[:,None]).tocsr();x.data=np.log1p(x.data)
            values[lo:hi]=x@weights+offset
    assert np.isfinite(values).all()
    source=m[['cell_id','dataset','capture_id','specimen_id','subject_key','reported_patient_link_supported']].copy()
    for j,name in enumerate(feature_names):source[name]=values[:,j]
    source.to_csv(dest,index=False)
    save_json(LOG/'module_cell_activity_definition.json',dict(cells=CELLS,modules=modules,representatives=reps.gene.tolist(),
        module_activity='Mean of member-gene Z scores using frozen full-pool mean and population SD of log1p(CP10K)',
        representative_values='Unstandardized log1p(CP10K) with nCount_input denominator',
        selection_or_modules_refitted=False,source_tables='All gene memberships and population means/SD saved',independent_validation_pending=True))
    return source,feature_names


def main():
    source,features=calculate_cells()
    out=TAB/'dynamic_genes/method_path_sensitivity';out.mkdir(exist_ok=False)
    summaries=[];roster=[]
    keys=['dataset','capture_id','specimen_id','subject_key','reported_patient_link_supported','bin']
    for method in ['Monocle2','Slingshot']:
        for root in ['primary','alternative']:
            suffix='path_cells' if method=='Monocle2' else 'cells'
            paths=pd.read_csv(TAB/f'{method}_{root}_{suffix}.csv.gz')
            assert np.array_equal(paths.cell_id,source.cell_id)
            for column in [x for x in paths if x.startswith('PT_')]:
                number=column.split('_')[1];pt=paths[column].to_numpy();w=paths['W_'+number].to_numpy()
                valid=np.isfinite(pt)&(w>=.5)
                name=f'{method}_{root}_Path{number}'
                row=dict(method=method,root=root,path='Path'+number,qualified_cells=int(valid.sum()))
                if valid.sum()<30 or np.ptp(pt[valid])==0:
                    row.update(status='insufficient_cells_or_time_range');roster.append(row);continue
                d=source.loc[valid].copy();time=pt[valid]
                d['bin']=np.minimum(((time-time.min())/np.ptp(time)*12).astype(int),11)
                capture=d.groupby(keys,observed=True).agg(cells=('cell_id','size'),**{f:(f,'mean') for f in features}).reset_index()
                capture['eligible']=capture.cells.ge(10);capture.to_csv(out/(name+'_capture_bins.csv.gz'),index=False)
                eligible=capture[capture.eligible]
                specimen=eligible.groupby(['dataset','specimen_id','subject_key','reported_patient_link_supported','bin'])[features].mean().reset_index()
                specimen.to_csv(out/(name+'_specimen_bins.csv.gz'),index=False)
                labels=specimen[specimen.reported_patient_link_supported].groupby(['dataset','subject_key','bin'])[features].mean().reset_index()
                labels.to_csv(out/(name+'_reported_label_bins.csv.gz'),index=False)
                proxy=specimen[~specimen.reported_patient_link_supported]
                proxy.to_csv(out/(name+'_unsupported_proxy_bins.csv.gz'),index=False)
                for b in range(12):
                    block=labels[labels.bin.eq(b)];unknown=proxy[proxy.bin.eq(b)]
                    for feature in features:
                        v=block[feature];u=unknown[feature]
                        summaries.append(dict(method=method,root=root,path='Path'+number,bin=b,feature=feature,
                            supported_reported_labels=len(block),supported=len(block)>=3,
                            raw_mean=float(v.mean()) if len(v) else np.nan,
                            mean=float(v.mean()) if len(v)>=3 else np.nan,
                            q25=float(v.quantile(.25)) if len(v)>=3 else np.nan,q75=float(v.quantile(.75)) if len(v)>=3 else np.nan,
                            unsupported_proxy_specimens=len(unknown),proxy_mean=float(u.mean()) if len(u) else np.nan,
                            proxy_min=float(u.min()) if len(u) else np.nan,proxy_max=float(u.max()) if len(u) else np.nan))
                row.update(status='descriptive',captures=capture.capture_id.nunique(),specimens=specimen.specimen_id.nunique(),
                    reported_labels=labels[['dataset','subject_key']].drop_duplicates().shape[0],unsupported_proxy_specimens=proxy.specimen_id.nunique())
                roster.append(row)
                print('Module path sensitivity',name,flush=True)
    summary=pd.DataFrame(summaries);summary.to_csv(out/'all_module_representative_path_summary.csv.gz',index=False)
    pd.DataFrame(roster).to_csv(out/'all_method_root_path_roster.csv',index=False)
    figout=ROOT/'results/01_图形/module_sensitivity';figout.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':8,'pdf.fonttype':42});manifest=[]
    limits={f:(float(d['mean'].min()),float(d['mean'].max())) for f,d in summary.groupby('feature')}
    for method in ['Monocle2','Slingshot']:
        for root in ['primary','alternative']:
            block=summary[summary.method.eq(method)&summary.root.eq(root)]
            paths=sorted(block.path.unique(),key=lambda x:int(x[4:]))
            fig,axes=plt.subplots(4,2,figsize=(15,max(12,1.2*len(paths))),layout='constrained')
            for ax,feature in zip(axes.flat,features):
                matrix=block[block.feature.eq(feature)].pivot(index='path',columns='bin',values='mean').reindex(index=paths,columns=range(12))
                cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#D6D6D6');lo,hi=limits[feature]
                im=ax.imshow(matrix,aspect='auto',interpolation='nearest',extent=(0,1,len(paths),0),cmap=cmap,vmin=lo,vmax=hi)
                ax.set(title=feature,yticks=np.arange(len(paths))+.5,yticklabels=paths,xlabel='Within-path scaled pseudotime')
                fig.colorbar(im,ax=ax,label='Mean member-gene Z' if feature.startswith('Module') else 'log1p(CP10K)',shrink=.8)
            fig.suptitle(f'{method} | {root} | every returned path\nEqual capture -> specimen -> supported reported label means (>=3 labels/bin); gray = insufficient coverage.\nFixed primary modules/representatives; proxies saved separately; no path selection or refitting.',fontsize=11)
            name=f'50_{method}_{root}_module_and_gene_paths'
            for ext in ['png','pdf']:fig.savefig(figout/(name+'.'+ext),dpi=180,bbox_inches='tight')
            plt.close(fig);manifest.append(dict(name=name,png='01_图形/module_sensitivity/'+name+'.png',pdf='01_图形/module_sensitivity/'+name+'.pdf',paths=paths,
                sources=['02_表格/AllEpithelial/dynamic_genes/method_path_sensitivity/all_module_representative_path_summary.csv.gz'],visual_review='pending'))
    save_json(LOG/'module_method_figure_manifest.json',manifest)
    save_json(LOG/'module_method_sensitivity_summary.json',dict(completed=True,all_paths=roster,
        feature_color_limits_across_all_methods_roots=limits,independent_numeric_validation_pending=True))


if __name__=='__main__':main()
