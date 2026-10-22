# Exercise the actual frozen engine and new entry on generated software data.
suppressPackageStartupMessages({library(slingshot);library(jsonlite)})
root=normalizePath(commandArgs(TRUE)[1])
fixture=file.path(root,'.runtime/exact_projection_e2e_fixture')
stopifnot(!dir.exists(fixture));dir.create(fixture,recursive=TRUE)
input=file.path(fixture,'work/AllEpithelial');dir.create(input,recursive=TRUE)
code=file.path(fixture,'results/code');dir.create(code,recursive=TRUE)
for(name in c('exact_projection_acceleration.R','exact_projection_bvh.cpp'))
 stopifnot(file.symlink(file.path(root,'results/code',name),file.path(code,name)))
stopifnot(file.symlink(file.path(root,'work/fit_generated.R'),file.path(fixture,'work/fit_generated.R')))
set.seed(20261020)
labels=rep(paste0('K',0:5),each=400)
centers=rbind(c(-4,0),c(-2,0),c(0,0),c(2,2),c(2,-2),c(4,-2))
xy=centers[match(labels,paste0('K',0:5)),]+matrix(rnorm(length(labels)*2,sd=.3),ncol=2)
basis=qr.Q(qr(matrix(rnorm(30*2),30,2)))
z=xy%*%t(basis)+matrix(rnorm(length(labels)*30,sd=.03),ncol=30)
colnames(z)=as.character(0:29)
metadata=data.frame(cell_id=paste0('synthetic_e2e_',seq_len(nrow(z))),coarse_state=labels,
 state_UMAP_1=xy[,1],state_UMAP_2=xy[,2])
write.csv(metadata,gzfile(file.path(input,'metadata.csv.gz')),row.names=FALSE)
write.csv(z,file.path(input,'Harmony.csv'),row.names=FALSE)
cfg=fromJSON(file.path(root,'config.json'));cfg$fitted_cells=nrow(z)
cfg$scope_name='Generated software fixtures, not biological cells'
cfg$blocks=data.frame(block='AllEpithelial',cells=nrow(z),root='K0',alternative='K5')
config=file.path(fixture,'work/analysis_request.json')
write_json(cfg,config,pretty=TRUE,auto_unbox=TRUE)
commands=list(native=file.path(root,'work/fit_generated.R'),
 indexed=file.path(root,'results/code/fit_slingshot_exact_entry.R'))
elapsed=list()
for(name in names(commands)){
 folder=file.path(fixture,name)
 for(child in c('objects','tables','logs'))dir.create(file.path(folder,child),recursive=TRUE)
 write_json(list(primary=list(actual_root_group='K0'),alternative=list(actual_root_group='K5')),file.path(folder,'tables/Monocle2_roots.json'),pretty=TRUE,auto_unbox=TRUE)
 args=c(commands[[name]],input,file.path(folder,'objects'),file.path(folder,'tables'),
        file.path(folder,'logs'),config,'Slingshot')
 started=proc.time()[['elapsed']]
 status=system2(Sys.which('Rscript'),args,
  stdout=file.path(folder,'console.log'),stderr=file.path(folder,'console.log'),
  env=c('OMP_NUM_THREADS=2','OPENBLAS_NUM_THREADS=2','MKL_NUM_THREADS=2','VECLIB_MAXIMUM_THREADS=2'))
 stopifnot(status==0);elapsed[[name]]=proc.time()[['elapsed']]-started
}
compare=function(a,b){
 stopifnot(identical(unname(is.na(a)),unname(is.na(b))),identical(dim(a),dim(b)))
 valid=is.finite(a)&is.finite(b);error=if(any(valid))max(abs(a[valid]-b[valid])) else 0
 stopifnot(error<1e-8);error
}
checks=list()
for(orientation in c('primary','alternative')){
 a=readRDS(file.path(fixture,'native/objects',paste0('Slingshot_',orientation,'.rds')))
 b=readRDS(file.path(fixture,'indexed/objects',paste0('Slingshot_',orientation,'.rds')))
 au=readRDS(file.path(fixture,'native/objects',paste0('Slingshot_',orientation,'_UMAP.rds')))
 bu=readRDS(file.path(fixture,'indexed/objects',paste0('Slingshot_',orientation,'_UMAP.rds')))
 stopifnot(identical(slingLineages(a),slingLineages(b)),
           identical(rownames(slingPseudotime(a)),metadata$cell_id),
           identical(rownames(slingPseudotime(b)),metadata$cell_id))
 checks[[orientation]]=list(cells=nrow(metadata),paths=ncol(slingPseudotime(a)),
   pseudotime_error=compare(slingPseudotime(a),slingPseudotime(b)),
   weights_error=compare(slingCurveWeights(a),slingCurveWeights(b)),
   projected_pseudotime_error=compare(slingPseudotime(au),slingPseudotime(bu)),
   curves_UMAP_error=max(vapply(seq_along(slingCurves(au)),function(i)
     compare(slingCurves(au)[[i]]$s,slingCurves(bu)[[i]]$s),numeric(1))))
}
audit=fromJSON(file.path(fixture,'indexed/logs/Slingshot_exact_projection_audit.json'))
stopifnot(nrow(audit)>=4L,all(audit$input_rows_retained==2400L))
report=list(passed=TRUE,scope='2400 generated software cells, never a biological subset',
 original_frozen_engine_used=TRUE,same_native_parameters=TRUE,
 comparison=checks,elapsed_seconds=elapsed,indexed_projection_calls=nrow(audit),
 formal_biological_data_used=FALSE,formal_biological_fit_restarted=FALSE)
write_json(report,file.path(root,'results/04_日志和检查/exact_entry_e2e_validation.json'),
 pretty=TRUE,auto_unbox=TRUE,digits=NA)
print(report)
