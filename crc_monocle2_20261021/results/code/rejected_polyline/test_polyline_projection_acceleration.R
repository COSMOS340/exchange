# Software fixtures only, never a subset of the CRC data.
suppressPackageStartupMessages({library(slingshot);library(jsonlite)})
root=normalizePath(commandArgs(TRUE)[1]);source(file.path(root,'results/code/rejected_polyline/polyline_projection_acceleration.R'))
accelerator=make_polyline_projection_accelerator(min_vertices=8L,tolerance=1e-12,max_vertices=256L)
native=accelerator$reference
compare=function(a,b,tolerance=1e-8){
 stopifnot(identical(unname(is.na(a)),unname(is.na(b))))
 valid=is.finite(a)&is.finite(b)
 stopifnot(all(is.finite(a)==is.finite(b)))
 error=if(any(valid))max(abs(a[valid]-b[valid])) else 0
 cat('comparison maximum error',format(error,digits=17),'\n')
 stopifnot(error<tolerance);error
}
set.seed(20261020)
piecewise=function(vertices,per_segment){
 do.call(rbind,lapply(seq_len(nrow(vertices)-1L),function(i){
  t=seq(0,1,length.out=per_segment)
  outer(1-t,vertices[i,])+outer(t,vertices[i+1,])
 }))
}
cases=list(
 straight_2D=piecewise(rbind(c(0,0),c(3,5)),250),
 branched_polyline_30D=piecewise(matrix(rnorm(4*30),4,30),80),
 backtracking=piecewise(rbind(c(0,0),c(2,1),c(0,0),c(4,-1)),80),
 duplicates=piecewise(rbind(c(0,0),c(0,0),c(1,2),c(1,2),c(3,1)),50),
 nearly_collinear_kink=piecewise(rbind(c(0,0),c(1,1e-9),c(2,0)),150))
records=list()
for(name in names(cases)){
 cat('FIXTURE',name,'\n')
 s=cases[[name]];x=matrix(rnorm(120*ncol(s)),120,ncol(s))
 a=native(x,s,stretch=0);b=accelerator$project(x,s,stretch=0)
 records[[name]]=list(projected_cells=nrow(x),vertices=nrow(s),
  coordinates_error=compare(a$s,b$s),pseudotime_error=compare(a$lambda,b$lambda),
  distance_error=compare(a$dist_ind,b$dist_ind))
}
# Verify that a slow accumulation of curvature cannot pass the global bound.
t=seq(0,1,length.out=2000)
curved=cbind(t,1e-7*t*t)
checked=accelerator$prune(curved,1e-12,256L)
stopifnot(!isTRUE(checked$safe))
stopifnot(!isTRUE(accelerator$prune(matrix(0,200,30),1e-12,256L)$safe))
# Same full synthetic cells, native settings, two computational reference roots.
set.seed(20261020)
labels=rep(paste0('K',0:5),each=60)
centers=rbind(c(-4,0),c(-2,0),c(0,0),c(2,2),c(2,-2),c(4,-2))
xy=centers[match(labels,paste0('K',0:5)),]+matrix(rnorm(length(labels)*2,sd=.3),ncol=2)
basis=qr.Q(qr(matrix(rnorm(30*2),30,2)))
z=xy%*%t(basis)+matrix(rnorm(length(labels)*30,sd=.03),ncol=30)
rownames(z)=paste0('synthetic_',seq_len(nrow(z)))
fits=list()
for(reference in c('K0','K5')){
 set.seed(20261020);s=getLineages(z,clusterLabels=labels,start.clus=reference)
 fits[[reference]]=getCurves(s,maxit=30,approx_points=150)
}
binding=accelerator$install()
comparisons=list()
for(reference in c('K0','K5')){
 set.seed(20261020);s=getLineages(z,clusterLabels=labels,start.clus=reference)
 s=getCurves(s,maxit=30,approx_points=150);baseline=fits[[reference]]
 stopifnot(identical(slingLineages(s),slingLineages(baseline)))
 comparisons[[reference]]=list(cells=nrow(z),dimensions=ncol(z),paths=ncol(slingPseudotime(s)),
  pseudotime_error=compare(slingPseudotime(s),slingPseudotime(baseline)),
  weight_error=compare(slingCurveWeights(s),slingCurveWeights(baseline)))
}
accelerator$uninstall()
stopifnot(length(accelerator$audit())>=length(cases))
report=list(passed=TRUE,scope='Synthetic software tests only; no CRC cells fitted or sampled',
 projection_cases=records,global_curvature_check_rejects=TRUE,degenerate_constant_curve_falls_back=TRUE,
 full_joint_synthetic_Slingshot=comparisons,patched_binding=binding,
 audit=accelerator$audit(),formal_biological_acceleration_applied=FALSE)
write_json(report,file.path(root,'results/04_日志和检查/polyline_acceleration_synthetic_validation.json'),
 pretty=TRUE,auto_unbox=TRUE,digits=NA)
print(report)
