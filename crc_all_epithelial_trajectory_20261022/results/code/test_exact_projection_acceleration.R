# Synthetic software verification against installed princurve/Slingshot.
suppressPackageStartupMessages({library(slingshot);library(jsonlite)})
root=normalizePath(commandArgs(TRUE)[1]);source(file.path(root,'results/code/exact_projection_acceleration.R'))
accelerator=make_exact_projection_accelerator(root,min_vertices=8L)
native=accelerator$reference
compare=function(a,b,tolerance=1e-8){
 stopifnot(identical(unname(is.na(a)),unname(is.na(b))),
           identical(unname(is.finite(a)),unname(is.finite(b))))
 valid=is.finite(a)&is.finite(b)
 error=if(any(valid))max(abs(a[valid]-b[valid])) else 0
 cat('maximum error',format(error,digits=17),'\n')
 stopifnot(error<tolerance);error
}
piecewise=function(vertices,per_segment){
 do.call(rbind,lapply(seq_len(nrow(vertices)-1L),function(i){
  t=seq(0,1,length.out=per_segment)
  outer(1-t,vertices[i,])+outer(t,vertices[i+1,])
 }))
}
set.seed(20261020)
cases=list(straight_2D=piecewise(rbind(c(0,0),c(3,5)),250),
 polyline_30D=piecewise(matrix(rnorm(4*30),4,30),80),
 backtracking=piecewise(rbind(c(0,0),c(2,1),c(0,0),c(4,-1)),80),
 duplicates=piecewise(rbind(c(0,0),c(0,0),c(1,2),c(1,2),c(3,1)),50),
 tiny_kink=piecewise(rbind(c(0,0),c(1,1e-9),c(2,0)),150),
 crossing=piecewise(rbind(c(-1,-1),c(1,1),c(-1,1),c(1,-1)),80),
 equidistant_arms=piecewise(rbind(c(-1,-1),c(-1,1),c(1,1),c(1,-1)),80))
for(i in 1:20){
 p=c(2,3,30)[(i-1)%%3+1]
 cases[[paste0('random_',i)]]=matrix(rnorm((30+i)*p),30+i,p)
}
records=list()
for(name in names(cases)){
 cat('FIXTURE',name,'\n')
 s=cases[[name]];x=matrix(rnorm(160*ncol(s)),160,ncol(s))
 if(name %in% c('crossing','equidistant_arms'))x=rbind(x,c(0,0),c(0,.5),c(0,-.5))
 rownames(x)=paste0('software_',seq_len(nrow(x)))
 a=native(x,s,stretch=0);b=accelerator$project(x,s,stretch=0)
 records[[name]]=list(rows=nrow(x),vertices=nrow(s),
  coordinates_error=compare(a$s,b$s),pseudotime_error=compare(a$lambda,b$lambda),
  distance_error=compare(a$dist_ind,b$dist_ind),ordering_identical=identical(a$ord,b$ord))
 stopifnot(identical(attributes(a),attributes(b)))
}
# Joint principal curves use all synthetic rows, two roots and native settings.
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
binding=accelerator$install();comparisons=list()
for(reference in c('K0','K5')){
 cat('JOINT_SYNTHETIC_ROOT',reference,'\n')
 set.seed(20261020);s=getLineages(z,clusterLabels=labels,start.clus=reference)
 s=getCurves(s,maxit=30,approx_points=150);baseline=fits[[reference]]
 stopifnot(identical(slingLineages(s),slingLineages(baseline)))
 comparisons[[reference]]=list(cells=nrow(z),dimensions=ncol(z),
  paths=ncol(slingPseudotime(s)),pseudotime_error=compare(slingPseudotime(s),slingPseudotime(baseline)),
  weight_error=compare(slingCurveWeights(s),slingCurveWeights(baseline)))
}
accelerator$uninstall()
# Timed original-vs-indexed projection on purely synthetic input.
set.seed(20261020)
s=piecewise(matrix(rnorm(4*30),4,30),2000)
x=matrix(rnorm(1000*30),1000,30)
original_time=system.time(a<-native(x,s,stretch=0))[['elapsed']]
indexed_time=system.time(b<-accelerator$project(x,s,stretch=0))[['elapsed']]
benchmark=list(rows=nrow(x),vertices=nrow(s),native_seconds=original_time,indexed_seconds=indexed_time,
 pseudotime_error=compare(a$lambda,b$lambda),coordinates_error=compare(a$s,b$s),
 distance_error=compare(a$dist_ind,b$dist_ind))
stopifnot(length(accelerator$audit())>length(cases))
report=list(passed=TRUE,scope='Synthetic software fixtures only; no CRC data sampled or fitted',
 projection_cases=records,joint_synthetic_Slingshot=comparisons,patched_binding=binding,
 benchmark=benchmark,audit=accelerator$audit(),all_original_segments_and_X_rows_retained=TRUE,
 formal_biological_acceleration_applied=FALSE)
write_json(report,file.path(root,'results/04_日志和检查/exact_projection_synthetic_validation.json'),
 pretty=TRUE,auto_unbox=TRUE,digits=NA)
print(report)
