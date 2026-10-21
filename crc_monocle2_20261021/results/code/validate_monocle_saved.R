# Independent validation of saved objects; never fits or changes the cell set.
suppressPackageStartupMessages({library(monocle);library(igraph);library(data.table);library(jsonlite)})
root=normalizePath(commandArgs(TRUE)[1],mustWork=TRUE)
tab=file.path(root,'results/02_表格/CRC');logs=file.path(root,'results/04_日志和检查')
metadata=as.data.frame(fread(file.path(root,'work/CRC/metadata.csv.gz')))
stopifnot(nrow(metadata)==150917L,!anyDuplicated(metadata$cell_id))

# Traverse exported tree edges directly, without calling igraph::distances.
edge_distances=function(g,root_id){
  ids=V(g)$name;n=length(ids);edges=as_edgelist(g,names=FALSE);w=E(g)$weight
  stopifnot(nrow(edges)==n-1L,all(is.finite(w)),all(w>=0))
  from=c(edges[,1],edges[,2]);to=c(edges[,2],edges[,1]);weight=rep(w,2)
  adj=split(seq_along(from),factor(from,levels=seq_len(n)))
  result=rep(NA_real_,n);start=match(root_id,ids);stopifnot(!is.na(start))
  result[start]=0;stack=integer(n);stack[1]=start;top=1L
  while(top>0L){
    u=stack[top];top=top-1L
    for(e in adj[[u]]){
      v=to[e]
      if(is.na(result[v])){
        result[v]=result[u]+weight[e];top=top+1L;stack[top]=v
      }
    }
  }
  stopifnot(all(is.finite(result)))
  setNames(result,ids)
}

canonical_edges=function(g){
  e=as_data_frame(g,what='edges')
  out=data.frame(a=pmin(e$from,e$to),b=pmax(e$from,e$to),weight=e$weight)
  out=out[order(out$a,out$b),];rownames(out)=NULL;out
}

checks=list();baseline=NULL
for(orientation in c('primary','alternative')){
  object=readRDS(file.path(root,'results/03_对象/CRC',paste0('Monocle2_',orientation,'.rds')))
  pd=pData(object);ids=as.character(pd$cell_id)
  stopifnot(identical(ids,metadata$cell_id),!anyDuplicated(ids),
            identical(as.character(pd$CNV_relative_band),metadata$CNV_relative_band),
            identical(as.numeric(pd$cnv_abs_log2),metadata$cnv_abs_log2))
  xy=reducedDimS(object);centers=reducedDimK(object);center_tree=minSpanningTree(object)
  stopifnot(ncol(xy)==150917L,all(is.finite(xy)),all(is.finite(centers)),
            is_tree(center_tree),is_connected(center_tree),
            ecount(center_tree)==vcount(center_tree)-1L)
  graph=object@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_tree
  actual_root=object@auxOrderingData[['DDRTree']]$root_cell
  stopifnot(is_tree(graph),is_connected(graph),vcount(graph)==150917L,
            setequal(V(graph)$name,ids),ecount(graph)==150916L)
  independent=edge_distances(graph,actual_root)[ids]
  error=max(abs(independent-as.numeric(pd$Pseudotime)))
  stopifnot(is.finite(error),error<1e-8)
  exported=fread(file.path(tab,paste0('Monocle2_',orientation,'_cells.csv.gz')))
  stopifnot(identical(exported$cell_id,ids),
            max(abs(exported$PT_1-independent))<1e-8,
            max(abs(as.matrix(exported[,.(DDRTree_1,DDRTree_2)])-t(xy)))<1e-8)
  geometry=list(xy=xy,centers=centers,edges=canonical_edges(center_tree))
  if(is.null(baseline))baseline=geometry else stopifnot(identical(baseline,geometry))
  fwrite(as_data_frame(graph,what='edges'),file.path(tab,paste0('Monocle2_',orientation,'_projected_tree_edges.csv.gz')))
  checks[[orientation]]=list(cells=length(ids),CNV_values_unchanged=TRUE,
    actual_root_cell=actual_root,actual_root_cluster=as.character(pd[actual_root,'coarse_state']),
    states=length(unique(pd$State)),center_vertices=vcount(center_tree),
    projected_vertices=vcount(graph),projected_edges=ecount(graph),
    independent_tree_distance_max_error=error,geometry_shared=TRUE)
  rm(object,pd,graph,geometry,xy,centers,center_tree,exported);gc()
}
objective=fread(file.path(tab,'Monocle2_objective.csv'))
checks$objective=list(entries=nrow(objective),finite=all(is.finite(objective$objective)),
  increasing_steps=sum(diff(objective$objective)>1e-8),
  last_absolute_change=if(nrow(objective)>1)tail(abs(diff(objective$objective)),1) else NA_real_)
checks$same_actual_root=identical(checks$primary$actual_root_cell,checks$alternative$actual_root_cell)
checks$passed=TRUE
checks$remaining='Figures, dynamics, method comparison, warning review and visual acceptance remain separate requirements'
write_json(checks,file.path(logs,'Monocle2_independent_validation.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
writeLines(capture.output(sessionInfo()),file.path(logs,'Monocle2_independent_validation_session.txt'))
print(checks)
