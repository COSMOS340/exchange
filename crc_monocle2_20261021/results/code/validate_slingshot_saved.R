# Reload both full-cell models and their UMAP projections without refitting.
suppressPackageStartupMessages({library(slingshot);library(igraph);library(data.table);library(jsonlite)})
root=normalizePath(commandArgs(TRUE)[1]);tab=file.path(root,'results/02_表格/CRC');log=file.path(root,'results/04_日志和检查')
m=as.data.frame(fread(file.path(root,'work/CRC/metadata.csv.gz')))
h=as.matrix(fread(file.path(root,'work/CRC/Harmony.csv')));rownames(h)=m$cell_id
cfg=fromJSON(file.path(root,'work/analysis_request.json'))
checks=list()
for(orientation in c('primary','alternative')){
  s=readRDS(file.path(root,'results/03_对象/CRC',paste0('Slingshot_',orientation,'.rds')))
  u=readRDS(file.path(root,'results/03_对象/CRC',paste0('Slingshot_',orientation,'_UMAP.rds')))
  pt=slingPseudotime(s);weights=slingCurveWeights(s)
  stopifnot(nrow(pt)==150917L,identical(rownames(pt),m$cell_id),all(is.finite(weights)),
    all(weights>=0),all(weights<=1+1e-12),identical(pt,slingPseudotime(u)),
    identical(weights,slingCurveWeights(u)))
  z=slingReducedDim(s)
  stopifnot(ncol(z)==30L,identical(rownames(z),m$cell_id),max(abs(z-h))<1e-10)
  d=fread(file.path(tab,paste0('Slingshot_',orientation,'_cells.csv.gz')))
  stopifnot(identical(d$cell_id,m$cell_id))
  errors=numeric(ncol(pt))
  for(j in seq_len(ncol(pt))){
    p=d[[paste0('PT_',j)]];w=d[[paste0('W_',j)]]
    stopifnot(identical(unname(is.na(p)),unname(is.na(pt[,j]))),max(abs(w-weights[,j]))<1e-8)
    finite=is.finite(pt[,j]);stopifnot(all(is.finite(p[finite])))
    errors[j]=if(any(finite))max(abs(p[finite]-pt[finite,j])) else 0
    stopifnot(errors[j]<1e-8)
  }
  g=slingMST(s);stopifnot(is_tree(g),is_connected(g))
  curves=slingCurves(u);saved=fread(file.path(tab,paste0('Slingshot_',orientation,'_curves.csv')))
  curve_error=0
  for(j in seq_along(curves)){
    expected=curves[[j]]$s[curves[[j]]$ord,1:2,drop=FALSE]
    actual=as.matrix(saved[saved$path==j,c('UMAP_1','UMAP_2'),with=FALSE])
    stopifnot(identical(dim(expected),dim(actual)))
    curve_error=max(curve_error,max(abs(expected-actual)))
  }
  stopifnot(curve_error<1e-8)
  requested=if(orientation=='primary')cfg$blocks$root[1] else cfg$blocks$alternative[1]
  lineages=slingLineages(s);stopifnot(all(vapply(lineages,function(x)x[1]==requested,logical(1))))
  checks[[orientation]]=list(cells=nrow(pt),dimensions=ncol(z),requested_root=requested,
    paths=ncol(pt),assigned_cells=sum(rowSums(is.finite(pt))>0),unassigned_cells=sum(rowSums(is.finite(pt))==0),
    pseudotime_export_max_error=max(errors),UMAP_curve_export_max_error=curve_error,
    UMAP_projection_preserves_pseudotime=TRUE,UMAP_projection_preserves_weights=TRUE,
    exact_cell_alignment=TRUE,all_original_weights_retained=TRUE)
  rm(s,u,pt,weights,z,d,g,curves,saved);gc()
}
checks$passed=TRUE
write_json(checks,file.path(log,'Slingshot_independent_validation.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
writeLines(capture.output(sessionInfo()),file.path(log,'Slingshot_final_session.txt'))
print(checks)
