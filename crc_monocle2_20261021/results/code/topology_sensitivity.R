# MST-only sensitivity. Inherited PCA/Harmony stay fixed; no curve/Harmony refit.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(slingshot);library(igraph)})
root=normalizePath(commandArgs(TRUE)[1]);input=file.path(root,'work/CRC');out=file.path(root,'results/02_表格/CRC')
m=as.data.frame(fread(file.path(input,'metadata.csv.gz')))
h=as.matrix(fread(file.path(input,'Harmony.csv')));p=as.matrix(fread(file.path(input,'PCA.csv')))
stopifnot(nrow(m)==150917L,nrow(h)==nrow(m),nrow(p)==nrow(m),ncol(h)==30L,ncol(p)==30L)
rownames(h)=rownames(p)=m$cell_id
request=fromJSON(file.path(root,'work/analysis_request.json'));reference=request$blocks$root[1]
edges=function(x,labels){
  s=getLineages(x,clusterLabels=labels,start.clus=if(reference %in% labels)reference else NULL)
  e=igraph::as_data_frame(slingMST(s))
  if(nrow(e)==0)return(character())
  sort(apply(e[,c('from','to'),drop=FALSE],1,function(z)paste(sort(z),collapse='|')))
}
main=edges(h,m$coarse_state)
records=list()
compare=function(label,e,n,clusters){
  union_n=length(union(main,e))
  data.frame(comparison=label,cells=n,clusters=clusters,main_edges=length(main),comparison_edges=length(e),
    shared_edges=length(intersect(main,e)),edge_Jaccard=if(union_n)length(intersect(main,e))/union_n else NA_real_,
    detail='MST-only sensitivity; inherited representations; no Harmony or curve refit')
}
records[[1]]=compare('PCA_vs_Harmony',edges(p,m$coarse_state),nrow(m),length(unique(m$coarse_state)))
for(ds in sort(unique(m$dataset))){
  keep=m$dataset!=ds
  e=tryCatch(edges(h[keep,,drop=FALSE],m$coarse_state[keep]),error=function(err)err)
  if(inherits(e,'error')){
    records[[length(records)+1]]=data.frame(comparison=paste0('leave_out_',ds),cells=sum(keep),clusters=length(unique(m$coarse_state[keep])),
      main_edges=length(main),comparison_edges=NA,shared_edges=NA,edge_Jaccard=NA,detail=conditionMessage(e))
  }else records[[length(records)+1]]=compare(paste0('leave_out_',ds),e,sum(keep),length(unique(m$coarse_state[keep])))
}
fwrite(rbindlist(records),file.path(out,'topology_sensitivity.csv'))
write_json(list(scope='MST-only, not refitted curves or Harmony',baseline_edges=main),file.path(out,'topology_sensitivity_scope.json'),pretty=TRUE,auto_unbox=TRUE)
