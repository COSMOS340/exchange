# Small cluster-MST diagnostic only: never fits or alters principal curves.
suppressPackageStartupMessages({library(slingshot);library(data.table);library(jsonlite)})
root=normalizePath(commandArgs(TRUE)[1])
m=as.data.frame(fread(file.path(root,'work/CRC/metadata.csv.gz')))
z=as.matrix(fread(file.path(root,'work/CRC/Harmony.csv'),header=TRUE))
stopifnot(nrow(z)==150917L,nrow(m)==150917L,ncol(z)==30L)
rownames(z)=m$cell_id
cfg=fromJSON(file.path(root,'work/analysis_request.json'))
counts=table(m$coarse_state)
records=list()
for(orientation in c('primary','alternative')){
 reference=if(orientation=='primary')cfg$blocks$root[1] else cfg$blocks$alternative[1]
 set.seed(20261020)
 s=getLineages(z,clusterLabels=m$coarse_state,start.clus=reference)
 paths=slingLineages(s)
 records[[orientation]]=list(reference=reference,lineages=paths,
   initial_path_cluster_member_cells=vapply(paths,function(x)sum(counts[x]),numeric(1)),
   initial_full_projection_pair_count=150917*sum(vapply(paths,function(x)sum(counts[x]),numeric(1))))
 rm(s);gc()
}
report=list(scope='MST initialization diagnostic only; no principal curves, no Harmony refit, no fitted model replacement',
 cells=nrow(z),dimensions=ncol(z),diagnostics=records,
 runtime_inference='Native getCurves performs full-cell projection onto all initial member projections before reducing to approx_points; pair counts describe initialization workload, not measured completion percentage')
write_json(report,file.path(root,'results/04_日志和检查/Slingshot_initial_MST_only_diagnostics.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
print(report)
