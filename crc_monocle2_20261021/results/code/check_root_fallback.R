# Synthetic geometry only; no biological cells, DDRTree fit or subsampling.
args=commandArgs(TRUE);root=normalizePath(args[1]);engine=file.path(root,'work/fit_generated.R')
text=paste(readLines(engine,warn=FALSE),collapse='\n')
finish=regexpr('options(warn=1)',text,fixed=TRUE)[1];stopifnot(finish>0)
eval(parse(text=substr(text,1,finish-1)))
patch_monocle2_dplyr();patch_monocle2_igraph();memory_patch();patch_root_fallback()
suppressPackageStartupMessages(library(jsonlite))
ids=paste0('v',1:15);genes=paste0('g',1:30)
counts=Matrix::Matrix(matrix(rep(1:5,length.out=30*15),30,15,
                            dimnames=list(genes,ids)),sparse=TRUE)
metadata=data.frame(cell_id=ids,State=factor(ifelse(1:15 %in% 7:9,'2','1')),
                    Pseudotime=seq(0,4,length.out=15),row.names=ids)
cds=newCellDataSet(as(counts,'dgCMatrix'),
  phenoData=new('AnnotatedDataFrame',data=metadata),
  featureData=new('AnnotatedDataFrame',data=data.frame(gene_short_name=genes,row.names=genes)),
  expressionFamily=negbinomial.size())
z=rbind(seq(-2,2,length.out=15),.02*sin(1:15));colnames(z)=ids
k=rbind(-2:2,rep(0,5));colnames(k)=paste0('Y_',1:5)
ns=asNamespace('monocle')
cds=get('reducedDimS<-',ns)(cds,value=z)
cds=get('reducedDimK<-',ns)(cds,value=k);cds@dim_reduce_type='DDRTree'
center_tree=graph_from_data_frame(data.frame(from=paste0('Y_',1:4),to=paste0('Y_',2:5),weight=1),directed=FALSE)
cds=get('minSpanningTree<-',ns)(cds,value=center_tree)
cds=get('cellPairwiseDistances<-',ns)(cds,value=as.matrix(dist(t(k))))
cds=get('findNearestPointOnMST',ns)(cds)
cds@auxOrderingData[['DDRTree']]$root_cell='v1'
selected_center=get('select_root_cell',asNamespace('monocle'))(cds,root_state='2')
stopifnot(selected_center %in% c('Y_2','Y_3','Y_4'))
messages=character()
# Use the patched namespace, matching the real pipeline's monocle::orderCells.
# The attached package binding can still point at the pre-patch function.
result=withCallingHandlers(monocle::orderCells(cds,root_state='2'),message=function(message){
  messages<<-c(messages,conditionMessage(message));invokeRestart('muffleMessage')})
stopifnot(any(grepl('ROOT_FALLBACK',messages,fixed=TRUE)))
graph=result@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_tree
leaves=V(graph)$name[degree(graph)==1]
center_index=match(selected_center,colnames(k))
expected_root=leaves[which.min(colSums((z[,leaves,drop=FALSE]-k[,center_index])^2))]
actual_root=result@auxOrderingData[['DDRTree']]$root_cell
stopifnot(identical(expected_root,actual_root),is_tree(graph),vcount(graph)==15L,
          identical(reducedDimS(result),z),identical(reducedDimK(result),k))
distance=as.numeric(distances(graph,v=actual_root,weights=E(graph)$weight))
error=max(abs(distance-pData(result)$Pseudotime));stopifnot(error<1e-8)
write_json(list(passed=TRUE,fixture='synthetic counts and a five-center line; internal requested root has no mapped projected leaf',
  synthetic_cells=15,selected_center=selected_center,expected_root=expected_root,actual_root=actual_root,
  pseudotime_distance_max_error=error,messages=messages,biological_fit_performed=FALSE),
  file.path(root,'results/04_日志和检查/root_fallback_synthetic_validation.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
