# Independent validation of saved objects; never fits or changes the cell set.
suppressPackageStartupMessages({library(monocle);library(igraph);library(data.table);library(jsonlite)})
root=normalizePath(commandArgs(TRUE)[1],mustWork=TRUE)
tab=file.path(root,'results/02_表格/CRC');logs=file.path(root,'results/04_日志和检查')
metadata=as.data.frame(fread(file.path(root,'work/CRC/metadata.csv.gz')))
stopifnot(nrow(metadata)==150917L,!anyDuplicated(metadata$cell_id))

# Traverse exported tree edges directly, without calling igraph::distances.
edge_distances=function(g,root_id){
  ids=V(g)$name;n=length(ids)
  # igraph returns double IDs. Integer IDs avoid factor coercing 100000 to
  # "1e+05" while explicit integer levels contain "100000".
  edges=matrix(as.integer(as_edgelist(g,names=FALSE)),ncol=2);w=E(g)$weight
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
  expected_sf=metadata$common_gene_UMI/exp(mean(log(metadata$common_gene_UMI)))
  sf_error=max(abs(as.numeric(pd$Size_Factor)-expected_sf)/expected_sf)
  stopifnot(all(is.finite(pd$Size_Factor)),all(pd$Size_Factor>0),sf_error<1e-10,
            identical(as.numeric(Matrix::colSums(exprs(object))),as.numeric(metadata$common_gene_UMI)))
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
  memberships=fread(file.path(tab,paste0('Monocle2_',orientation,'_path_cells.csv.gz')))
  definition=fromJSON(file.path(tab,paste0('Monocle2_',orientation,'_paths.json')),simplifyVector=FALSE)
  nearest=object@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_closest_vertex
  center_ids=V(center_tree)$name
  stopifnot(identical(memberships$cell_id,ids))
  for(j in seq_along(definition$paths)){
    inside=center_ids[nearest[,1]] %in% unlist(definition$paths[[j]])
    stopifnot(identical(as.logical(memberships[[paste0('W_',j)]]),as.logical(inside)),
              identical(is.finite(memberships[[paste0('PT_',j)]]),as.logical(inside)),
              max(abs(memberships[[paste0('PT_',j)]][inside]-pd$Pseudotime[inside]))<1e-8)
  }
  geometry=list(xy=xy,centers=centers,edges=canonical_edges(center_tree))
  if(is.null(baseline))baseline=geometry else stopifnot(identical(baseline,geometry))
  fwrite(as_data_frame(graph,what='edges'),file.path(tab,paste0('Monocle2_',orientation,'_projected_tree_edges.csv.gz')))
  checks[[orientation]]=list(cells=length(ids),CNV_values_unchanged=TRUE,
    actual_root_cell=actual_root,actual_root_cluster=as.character(pd[actual_root,'coarse_state']),
    states=length(unique(pd$State)),center_vertices=vcount(center_tree),
    projected_vertices=vcount(graph),projected_edges=ecount(graph),
    independent_tree_distance_max_error=error,geometry_shared=TRUE,
    independent_size_factor_max_relative_error=sf_error,
    path_membership_independently_verified=TRUE,paths=length(definition$paths))
  if(orientation=='primary'){
    disp=dispersionTable(object)
    fwrite(disp,file.path(tab,'Monocle2_dispersion_estimates.csv'))
    checks$dispersion=list(genes=nrow(disp),
      fitted_nonfinite=sum(!is.finite(disp$dispersion_fit)),
      empirical_nonfinite=sum(!is.finite(disp$dispersion_empirical)),
      console_warnings_require_review=TRUE)
  }
  rm(object,pd,graph,geometry,xy,centers,center_tree,exported);gc()
}
objective=fread(file.path(tab,'Monocle2_objective.csv'))
checks$objective=list(entries=nrow(objective),finite=all(is.finite(objective$objective)),
  increasing_steps=sum(diff(objective$objective)>1e-8),
  last_absolute_change=if(nrow(objective)>1)tail(abs(diff(objective$objective)),1) else NA_real_,
  last_relative_change=if(nrow(objective)>1)abs(tail(objective$objective,1)-tail(objective$objective,2)[1])/abs(tail(objective$objective,2)[1]) else NA_real_)
stopifnot(checks$objective$finite)
checks$same_actual_root=identical(checks$primary$actual_root_cell,checks$alternative$actual_root_cell)
checks$passed=TRUE
checks$remaining='Figures, dynamics, method comparison, warning review and visual acceptance remain separate requirements'
write_json(checks,file.path(logs,'Monocle2_independent_validation.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA)
writeLines(capture.output(sessionInfo()),file.path(logs,'Monocle2_independent_validation_session.txt'))
print(checks)
