
patch_blocked_dispersion <- function(){
 ns=asNamespace('monocle');original=get('disp_calc_helper_NB',ns)
 fast=function(cds,expressionFamily,min_cells_detected){
   counts=Biobase::exprs(cds);sf=Biobase::pData(cds)$Size_Factor
   if(!inherits(counts,'dgCMatrix') || any(!is.finite(sf)) || any(sf<=0))
     return(original(cds,expressionFamily,min_cells_detected))
   sums=squares=detected=numeric(nrow(counts))
   for(start in seq.int(1L,ncol(counts),by=1024L)){
     ids=seq.int(start,min(ncol(counts),start+1023L))
     block=round(counts[,ids,drop=FALSE])
     detected=detected+Matrix::rowSums(block>cds@lowerDetectionLimit)
     block@x=block@x/rep.int(sf[ids],diff(block@p))
     sums=sums+Matrix::rowSums(block)
     squares=squares+Matrix::rowSums(block^2)
   }
   keep=which(detected>min_cells_detected)
   mu=as.numeric(sums[keep]/ncol(counts))
   variance=pmax(0,as.numeric(squares[keep]/ncol(counts))-mu^2)
   dispersion=(variance-mean(1/sf)*mu)/mu^2
   dispersion[dispersion<0]=0;dispersion[mu==0]=NA;mu[mu==0]=NA
   data.frame(gene_id=rownames(counts)[keep],mu=mu,disp=dispersion)
 }
 assignInNamespace('disp_calc_helper_NB',fast,ns='monocle')
 # Avoid copying the complete expression matrix for an unchanged column order.
 env=new.env(parent=ns)
 env$dispersion_subset=function(cds,ids){
   if(identical(as.character(ids),colnames(cds)))cds else cds[,ids]
 }
 replacements=0L
 rewrite=function(x){
   if(identical(x,quote(cds[, .$rowname]))){
     replacements<<-replacements+1L;return(quote(dispersion_subset(cds,.$rowname)))
   }
   if(is.call(x))for(i in seq_along(x))
     if(!identical(x[[i]],quote(expr=)))x[[i]]=rewrite(x[[i]])
   x
 }
 estimator=get('estimateDispersionsForCellDataSet',ns)
 body(estimator)=rewrite(body(estimator));stopifnot(replacements==2L)
 environment(estimator)=env
 assignInNamespace('estimateDispersionsForCellDataSet',estimator,ns='monocle')
}


patch_sparse_size_factors <- function() {
 ns=asNamespace('monocle'); original=get('estimateSizeFactorsForSparseMatrix',ns)
 stopifnot(grepl('sfs <- cell_total/exp(mean(log(cell_total)))',
                paste(deparse(original),collapse='\n'),fixed=TRUE))
 fast=function(counts,locfunc=median,round_exprs=TRUE,method='mean-geometric-mean-total'){
   if(method!='mean-geometric-mean-total' || !inherits(counts,'dgCMatrix'))
     return(original(counts,locfunc=locfunc,round_exprs=round_exprs,method=method))
   total=numeric(ncol(counts));names(total)=colnames(counts)
   for(start in seq.int(1L,ncol(counts),by=1024L)){
     ids=seq.int(start,min(ncol(counts),start+1023L))
     block=counts[,ids,drop=FALSE]
     if(round_exprs)block=round(block)
     total[ids]=Matrix::colSums(block)
   }
   sfs=total/exp(mean(log(total)));sfs[is.na(sfs)]=1
   sfs
 }
 assignInNamespace('estimateSizeFactorsForSparseMatrix',fast,ns='monocle')
}

patch_monocle2_dplyr <- function() {
    select_compat <- function(.data, ..., .dots = list()) {
        specs <- c(list(...), as.list(.dots))
        if (length(specs) == 0) 
            return(.data[, integer(0), drop = FALSE])
        selected <- character(length(specs))
        output_names <- names(specs)
        if (is.null(output_names)) 
            output_names <- rep("", length(specs))
        for (i in seq_along(specs)) {
            spec <- specs[[i]]
            if (is.numeric(spec)) {
                selected[i] <- names(.data)[as.integer(spec)[1]]
            }
            else {
                selected[i] <- as.character(spec)[1]
            }
            if (!nzchar(output_names[i])) 
                output_names[i] <- selected[i]
        }
        result <- .data[, selected, drop = FALSE]
        names(result) <- output_names
        result
    }
    group_by_compat <- function(.data, ..., .dots = list(), add = FALSE) {
        specs <- c(list(...), as.list(.dots))
        groups <- unlist(lapply(specs, function(x) as.character(x)[1]), 
            use.names = FALSE)
        groups <- groups[nzchar(groups)]
        if (length(groups) == 0) 
            return(dplyr::group_by(.data, .add = add))
        dplyr::group_by(.data, !!!rlang::syms(groups), .add = add)
    }
    environment(select_compat) <- asNamespace("dplyr")
    environment(group_by_compat) <- asNamespace("dplyr")
    assignInNamespace("select_", select_compat, ns = "dplyr")
    assignInNamespace("group_by_", group_by_compat, ns = "dplyr")
    monocle_imports <- parent.env(asNamespace("monocle"))
    replace_monocle_import <- function(name, replacement) {
        if (!exists(name, envir = monocle_imports, inherits = FALSE)) 
            return(invisible(FALSE))
        was_locked <- bindingIsLocked(name, monocle_imports)
        if (was_locked) 
            unlockBinding(name, monocle_imports)
        assign(name, replacement, envir = monocle_imports)
        if (was_locked) 
            lockBinding(name, monocle_imports)
        invisible(TRUE)
    }
    replace_monocle_import("select_", select_compat)
    replace_monocle_import("group_by_", group_by_compat)
}
patch_monocle2_igraph <- function() {
    ns <- asNamespace("monocle")
    project2mst_compat <- function(cds, Projection_Method) {
        dp_mst <- minSpanningTree(cds)
        Z <- reducedDimS(cds)
        Y <- reducedDimK(cds)
        cds <- findNearestPointOnMST(cds)
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex_names <- colnames(Y)[closest_vertex]
        closest_vertex_df <- as.matrix(closest_vertex)
        row.names(closest_vertex_df) <- row.names(closest_vertex)
        tip_leaves <- names(which(igraph::degree(dp_mst) == 1))
        if (!is.function(Projection_Method)) {
            P <- Y[, closest_vertex]
        }
        else {
            P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
            for (i in seq_along(closest_vertex)) {
                neighbors <- igraph::as_ids(igraph::neighbors(dp_mst, 
                  v = closest_vertex_names[i], mode = "all"))
                projection <- NULL
                distance <- NULL
                Z_i <- Z[, i]
                for (neighbor in neighbors) {
                  if (closest_vertex_names[i] %in% tip_leaves) {
                    tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], 
                      neighbor)])
                  }
                  else {
                    tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], 
                      neighbor)])
                  }
                  projection <- rbind(projection, tmp)
                  distance <- c(distance, stats::dist(rbind(Z_i, 
                    tmp)))
                }
                if (!methods::is(projection, "matrix")) {
                  projection <- as.matrix(projection)
                }
                P[, i] <- projection[which(distance == min(distance))[1], 
                  ]
            }
        }
        colnames(P) <- colnames(Z)
        dp <- as.matrix(stats::dist(t(P)))
        min_dist <- min(dp[dp != 0])
        dp <- dp + min_dist
        diag(dp) <- 0
        cellPairwiseDistances(cds) <- dp
        gp <- igraph::graph_from_adjacency_matrix(dp, mode = "undirected", 
            weighted = TRUE, diag = FALSE)
        dp_mst <- igraph::mst(gp)
        cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
        cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
        cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
        cds
    }
    environment(project2mst_compat) <- ns
    assignInNamespace("project2MST", project2mst_compat, ns = "monocle")
    extract_ordering_compat <- function(cds, root_cell, verbose = TRUE) {
        dp <- cellPairwiseDistances(cds)
        dp_mst <- minSpanningTree(cds)
        current_state <- 1
        states <- rep(1, ncol(dp))
        names(states) <- igraph::V(dp_mst)$name
        pseudotimes <- rep(0, ncol(dp))
        names(pseudotimes) <- igraph::V(dp_mst)$name
        parents <- rep(NA_character_, ncol(dp))
        names(parents) <- igraph::V(dp_mst)$name
        traversal <- igraph::dfs(dp_mst, root = root_cell, mode = "all", 
            unreachable = FALSE, parent = TRUE)
        traversal$parent <- as.numeric(traversal$parent)
        for (i in seq_along(traversal$order)) {
            current_node <- as.numeric(traversal$order[i])
            current_name <- igraph::V(dp_mst)[current_node]$name
            if (!is.na(traversal$parent[current_node])) {
                parent_node <- traversal$parent[current_node]
                parent_name <- igraph::V(dp_mst)[parent_node]$name
                current_pseudotime <- pseudotimes[parent_name] + 
                  dp[current_name, parent_name]
                if (igraph::degree(dp_mst, v = parent_name) > 
                  2) {
                  current_state <- current_state + 1
                }
            }
            else {
                parent_name <- NA_character_
                current_pseudotime <- 0
            }
            pseudotimes[current_name] <- current_pseudotime
            states[current_name] <- current_state
            parents[current_name] <- parent_name
        }
        result <- data.frame(sample_name = names(states), cell_state = factor(states), 
            pseudo_time = as.vector(pseudotimes), parent = parents)
        row.names(result) <- result$sample_name
        result
    }
    environment(extract_ordering_compat) <- ns
    assignInNamespace("extract_ddrtree_ordering", extract_ordering_compat, 
        ns = "monocle")
}


suppressPackageStartupMessages({library(igraph);library(Rcpp)})
Rcpp::cppFunction('List explicit_dfs(List adj, int root) {
 int n=adj.size(), count=0;
 IntegerVector parent(n,NA_INTEGER), order(n,NA_INTEGER), next(n,0);
 std::vector<bool> seen(n,false); std::vector<int> stack;
 stack.push_back(root-1); seen[root-1]=true; order[count++]=root;
 while(!stack.empty()){
   int u=stack.back(); IntegerVector neighbors=adj[u];
   if(next[u]>=neighbors.size()){stack.pop_back();continue;}
   int v=neighbors[next[u]++]-1;
   if(!seen[v]){seen[v]=true;parent[v]=u+1;order[count++]=v+1;stack.push_back(v);}
 }
 if(count!=n)stop("Trajectory tree is disconnected");
 return List::create(_["parent"]=parent,_["order"]=order);
}')
iterative_dfs=function(g,root){
 stopifnot(!is_directed(g),is_tree(g))
 index=if(is.character(root))match(root,V(g)$name) else as.integer(root)
 stopifnot(length(index)==1,!is.na(index))
 explicit_dfs(lapply(as_adj_list(g,mode="all"),as.integer),index)
}


patch_root_fallback=function(){
 ns=asNamespace("monocle");f=get("orderCells",ns)
 txt=paste(deparse(f),collapse="\n")
 old=paste0("        if (is.na(root_cell)) {\n",
            "            root_cell <- select_root_cell(cds, root_state, reverse)\n",
            "        }")
 stopifnot(grepl(old,txt,fixed=TRUE))
 replacement=paste0("        if (is.na(root_cell)) {\n",
   "            stopifnot(length(root_cell_idx)==1,length(tip_leaves)>0)\n",
   "            delta <- reducedDimS(cds)[, tip_leaves, drop=FALSE] - K_old[, root_cell_idx]\n",
   "            root_cell <- tip_leaves[which.min(colSums(delta^2))]\n",
   "            message('ROOT_FALLBACK: selected center ', V(old_mst)$name[root_cell_idx], ",
   "                      ' has no projected leaf; nearest DDRTree-space leaf is ', root_cell)\n",
   "        }")
 txt=sub(old,replacement,txt,fixed=TRUE)
 f=eval(parse(text=txt));environment(f)=ns
 assignInNamespace("orderCells",f,ns="monocle")
}


suppressPackageStartupMessages({library(monocle);library(igraph);library(Rcpp)})
Rcpp::cppFunction('List complete_prim(NumericMatrix z, bool skip_zero=false) {
 int n=z.nrow(), d=z.ncol();
 if(n<2) stop("At least two vertices required");
 std::vector<double> best(n, R_PosInf);
 std::vector<int> parent(n, -1); std::vector<bool> used(n,false);
 IntegerVector from(n-1), to(n-1); NumericVector weight(n-1);
 best[0]=0; double minpos=R_PosInf;
 for(int it=0;it<n;++it){
   int u=-1; double v=R_PosInf;
   for(int j=0;j<n;++j) if(!used[j] && best[j]<v){v=best[j];u=j;}
   if(u<0) stop("Disconnected complete graph under zero-edge policy");
   used[u]=true;
   if(it>0){from[it-1]=parent[u]+1;to[it-1]=u+1;weight[it-1]=std::sqrt(v);}
   for(int j=0;j<n;++j) if(!used[j]){
     double ds=0;for(int k=0;k<d;++k){double a=z(u,k)-z(j,k);ds+=a*a;}
     if(ds>0 && ds<minpos)minpos=ds;
     if((!skip_zero || ds>0) && ds<best[j]){best[j]=ds;parent[j]=u;}
   }
   if(it%500==0) Rcpp::checkUserInterrupt();
 }
 return List::create(_["from"]=from,_["to"]=to,_["weight"]=weight,
                     _["min_positive_distance"]=std::sqrt(minpos));
}')
complete_tree <- function(z, shifted=FALSE, skip_zero=FALSE) {
  stopifnot(all(is.finite(z)), !is.null(rownames(z)))
  e=complete_prim(z,skip_zero)
  if(shifted){stopifnot(is.finite(e$min_positive_distance));e$weight=e$weight+e$min_positive_distance}
  igraph::graph_from_data_frame(data.frame(from=rownames(z)[e$from],
      to=rownames(z)[e$to],weight=e$weight),directed=FALSE,
      vertices=data.frame(name=rownames(z)))
}
memory_patch <- function() {
 ns=asNamespace('monocle')
 patch_env=new.env(parent=ns);patch_env$complete_tree=complete_tree
 # E[(X-E[X])^2] = E[X^2]-E[X]^2, evaluated without densifying sparse RNA.
 f=get('disp_calc_helper_NB',ns);txt=paste(deparse(f),collapse='\n')
 old='f_expression_var <- Matrix::rowMeans((x - f_expression_mean)^2)'
 stopifnot(grepl(old,txt,fixed=TRUE))
 txt=sub(old,'f_expression_var <- pmax(0, Matrix::rowMeans(x^2) - as.numeric(f_expression_mean)^2)',txt,fixed=TRUE)
 f=eval(parse(text=txt));environment(f)=ns
 assignInNamespace('disp_calc_helper_NB',f,ns='monocle')
 # Start with the audited compatibility function, replacing only graph storage.
 f=get('project2MST',ns);txt=paste(deparse(f),collapse='\n')
 a=regexpr('    dp <- as.matrix',txt,fixed=TRUE)[1]
 b=regexpr('    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree',txt,fixed=TRUE)[1]
 stopifnot(a>0,b>a)
 txt=paste0(substr(txt,1,a-1),
   '    dp_mst <- complete_tree(t(P), shifted=TRUE)\n',
   '    cellPairwiseDistances(cds) <- matrix(numeric(0),0,0)\n',substr(txt,b,nchar(txt)))
 f=eval(parse(text=txt));environment(f)=patch_env
 assignInNamespace('project2MST',f,ns='monocle')
 extract=function(cds,root_cell,verbose=TRUE){
   g=minSpanningTree(cds);ids=igraph::V(g)$name;n=length(ids)
   tr=iterative_dfs(g,root_cell)
   par=as.integer(tr$parent);ord=as.integer(tr$order);deg=igraph::degree(g)
   pt=as.numeric(igraph::distances(g,v=root_cell,weights=igraph::E(g)$weight))
   state=integer(n);st=1L
   for(i in ord){if(!is.na(par[i]) && deg[par[i]]>2)st=st+1L;state[i]=st}
   data.frame(sample_name=ids,cell_state=factor(state),pseudo_time=pt,
     parent=ids[par],row.names=ids)
 }
 assignInNamespace('extract_ddrtree_ordering',extract,ns='monocle')
 f=get('select_root_cell',ns);txt=paste(deparse(f),collapse='\n')
 a=regexpr('        dp <- as.matrix',txt,fixed=TRUE)[1]
 b=regexpr('        tip_leaves <-',txt,fixed=TRUE)[1]
 stopifnot(a>0,b>a)
 txt=paste0(substr(txt,1,a-1),
   '        dp_mst <- complete_tree(t(reducedDimS(cds)[, row.names(root_cell_candidates),drop=FALSE]), skip_zero=TRUE)\n',
   substr(txt,b,nchar(txt)))
 f=eval(parse(text=txt));environment(f)=patch_env
 assignInNamespace('select_root_cell',f,ns='monocle')
}

options(warn=1)

suppressPackageStartupMessages({library(data.table);library(jsonlite);library(Matrix)})
args=commandArgs(TRUE);input=args[1];output=args[2];tab=args[3];logdir=args[4];config=args[5];method=args[6]
dir.create(output,recursive=TRUE,showWarnings=FALSE);dir.create(tab,recursive=TRUE,showWarnings=FALSE)
cfg=fromJSON(config);key=basename(input);spec=cfg$blocks[cfg$blocks$block==key,,drop=FALSE]
m=as.data.frame(fread(file.path(input,'metadata.csv.gz')));rownames(m)=m$cell_id
status=function(stage){
 cat(format(Sys.time()),key,method,stage,'cells',nrow(m),'\n');flush.console()
 write_json(list(block=key,method=method,stage=stage,cells=nrow(m),pid=Sys.getpid(),
   updated_at=format(Sys.time(),'%Y-%m-%d %H:%M:%S')),file.path(logdir,paste0(key,'_',method,'_progress.json')),pretty=TRUE,auto_unbox=TRUE)
}
status('loading');start=proc.time()[['elapsed']]
writeLines(capture.output(sessionInfo()),file.path(logdir,paste0(key,'_',method,'_session.txt')))
if(method=='Harmony'){
 suppressPackageStartupMessages(library(harmony))
 x=as.matrix(fread(file.path(input,'PCA.csv'),header=TRUE))
 stopifnot(nrow(x)==nrow(m),ncol(x)==30)
 h=RunHarmony(x,m,vars_use='dataset',theta=2,lambda=1,nclust=50,max_iter=20,ncores=2)
 fwrite(as.data.table(h),file.path(input,'Harmony.csv'))
 write_json(list(method=method,block=key,cells=nrow(m),all_cells_fitted=TRUE,sampled_cells=0),
   file.path(logdir,paste0(key,'_',method,'_completed.json')),pretty=TRUE,auto_unbox=TRUE)
 status('complete');quit(save='no')
}
if(method=='Slingshot'){
 suppressPackageStartupMessages(library(slingshot))
 z=as.matrix(fread(file.path(input,'Harmony.csv'),header=TRUE));rownames(z)=m$cell_id
 u=as.matrix(m[,c('state_UMAP_1','state_UMAP_2')]);rownames(u)=m$cell_id
 stopifnot(ncol(z)==30,nrow(z)==nrow(m),all(is.finite(z)),all(is.finite(u)))
 roots=list(primary=spec$root,alternative=spec$alternative)
 models=list();summary=list()
 for(name in names(roots)){
   root=roots[[name]];if(length(root)==0 || is.na(root))root=NULL
   if(name=='alternative' && length(unique(m$coarse_state))==1)next
   status(paste0('fit_',name))
   set.seed(20261020)
   s=getLineages(z,clusterLabels=m$coarse_state,start.clus=root)
   s=getCurves(s,maxit=cfg$slingshot_max_iter,approx_points=cfg$curve_points)
   pt=slingPseudotime(s);w=slingCurveWeights(s)
   stopifnot(nrow(pt)==nrow(m),identical(rownames(pt),m$cell_id))
   result=data.frame(cell_id=m$cell_id)
   for(i in seq_len(ncol(pt))){result[[paste0('PT_',i)]]=pt[,i];result[[paste0('W_',i)]]=w[,i]}
   fwrite(result,file.path(tab,paste0('Slingshot_',name,'_cells.csv.gz')))
   saveRDS(s,file.path(output,paste0('Slingshot_',name,'.rds')),compress=FALSE)
   before=pt;su=embedCurves(s,u)
   stopifnot(identical(before,slingPseudotime(su)))
   saveRDS(su,file.path(output,paste0('Slingshot_',name,'_UMAP.rds')),compress=FALSE)
   curves=slingCurves(su)
   ct=rbindlist(lapply(seq_along(curves),function(i){
      v=curves[[i]];data.frame(path=i,point=seq_along(v$ord),
         UMAP_1=v$s[v$ord,1],UMAP_2=v$s[v$ord,2])}))
   fwrite(ct,file.path(tab,paste0('Slingshot_',name,'_curves.csv')))
   ed=igraph::as_data_frame(slingMST(s));fwrite(ed,file.path(tab,paste0('Slingshot_',name,'_MST.csv')))
   summary[[name]]=list(requested_root=root,paths=slingLineages(s),assigned_cells=sum(rowSums(is.finite(pt))>0),
      unassigned_cells=sum(rowSums(is.finite(pt))==0),UMAP_projection_preserves_pseudotime=TRUE,
      path_cells=colSums(is.finite(pt)),dimensions=ncol(z),fit_all_cells=nrow(m))
 }
 if(length(unique(m$coarse_state))==1){
    result=as.data.frame(fread(file.path(tab,'Slingshot_primary_cells.csv.gz')))
    result$PT_1=max(result$PT_1,na.rm=TRUE)-result$PT_1
    fwrite(result,file.path(tab,'Slingshot_alternative_cells.csv.gz'))
    summary$alternative=list(type='orientation reversal only; same one-cluster curve, not an independently refitted root model')
 }
 write_json(summary,file.path(tab,'Slingshot_roots_and_paths.json'),pretty=TRUE,auto_unbox=TRUE)
}else if(method=='Monocle2'){
 patch_monocle2_dplyr();patch_monocle2_igraph();memory_patch();patch_root_fallback();patch_sparse_size_factors();patch_blocked_dispersion()
 checkpoint=file.path(output,'Monocle2_DDRTree_checkpoint.rds')
 if(file.exists(checkpoint)){
   status('resume_saved_DDRTree');cds=readRDS(checkpoint)
 }else{
   h=file.path(input,'counts.h5');shape=as.integer(rhdf5::h5read(h,'shape'))
   x=new('dgCMatrix',x=as.numeric(rhdf5::h5read(h,'data')),
      i=as.integer(rhdf5::h5read(h,'indices')),p=as.integer(rhdf5::h5read(h,'indptr')),Dim=shape)
   genes=as.data.frame(fread(file.path(input,'genes.csv')))
   rownames(x)=genes$gene;colnames(x)=m$cell_id
   stopifnot(ncol(x)==nrow(m),all(x@x>=0),all(x@x==floor(x@x)))
   cds=newCellDataSet(x,phenoData=new('AnnotatedDataFrame',data=m),
     featureData=new('AnnotatedDataFrame',data=data.frame(gene_short_name=genes$gene,row.names=genes$gene)),
     expressionFamily=negbinomial.size());rm(x);gc()
   status('size_factors');cds=estimateSizeFactors(cds);gc();status('dispersion');cds=estimateDispersions(cds,cores=1);gc()
   cds=detectGenes(cds,min_expr=0.1);cds=setOrderingFilter(cds,genes$gene[genes$ordering])
   stopifnot(all(fData(cds)$num_cells_expressed[genes$ordering]>=10))
   status('DDRTree_all_cells')
   cds=reduceDimension(cds,max_components=2,reduction_method='DDRTree',
      residualModelFormulaStr='~dataset',maxIter=cfg$monocle_max_iter,tol=cfg$monocle_tolerance)
   saveRDS(cds,checkpoint,compress=FALSE)
 }
 coord=t(reducedDimS(cds));colnames(coord)=c('DDRTree_1','DDRTree_2')
 centers=t(reducedDimK(cds));colnames(centers)=c('DDRTree_1','DDRTree_2')
 fwrite(data.frame(center=rownames(centers),centers),file.path(tab,'Monocle2_centers.csv'))
 fwrite(igraph::as_data_frame(minSpanningTree(cds)),file.path(tab,'Monocle2_tree_edges.csv'))
 objvals=cds@auxOrderingData$DDRTree$objective_vals
 fwrite(data.frame(iteration=seq_along(objvals),objective=as.numeric(objvals)),file.path(tab,'Monocle2_objective.csv'))
 status('complete_graph_ordering');base=monocle::orderCells(cds)
 choose=function(label,object){
   if(length(label)==0 || is.na(label))return(NULL)
   d=as.data.frame(table(pData(object)$State,pData(object)$coarse_state))
   colnames(d)=c('state','label','n');d=d[d$label==label,]
   total=table(pData(object)$State);d$fraction=d$n/as.numeric(total[as.character(d$state)])
   d=d[order(-d$fraction,-d$n,as.integer(as.character(d$state))),]
   as.character(d$state[1])
 }
 summaries=list()
 for(name in c('primary','alternative')){
    requested=if(name=='primary')spec$root else spec$alternative
    state=choose(requested,base)
    status(paste0('root_',name))
    if(is.null(state))res=monocle::orderCells(cds,reverse=name=='alternative') else res=monocle::orderCells(base,root_state=state)
    stopifnot(identical(reducedDimS(res),reducedDimS(cds)),identical(reducedDimK(res),reducedDimK(cds)))
    pd=pData(res);root=res@auxOrderingData$DDRTree$root_cell
    g=res@auxOrderingData$DDRTree$pr_graph_cell_proj_tree
    independent=as.numeric(igraph::distances(g,v=root,weights=igraph::E(g)$weight))
    stopifnot(all(is.finite(pd$Pseudotime)),max(abs(independent-pd$Pseudotime))<1e-8)
    fwrite(data.frame(cell_id=m$cell_id,coord,State=pd$State,PT_1=pd$Pseudotime,W_1=1),
         file.path(tab,paste0('Monocle2_',name,'_cells.csv.gz')))
    center_tree=minSpanningTree(res)
    nearest=res@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex
    center_ids=V(center_tree)$name
    center_root=center_ids[nearest[root,1]]
    leaves=setdiff(names(which(degree(center_tree)==1)),center_root)
    allpaths=igraph::shortest_paths(center_tree,from=center_root,to=leaves,output='vpath')$vpath
    memberships=data.frame(cell_id=m$cell_id);pdefs=list()
    for(j in seq_along(allpaths)){
       ids=as_ids(allpaths[[j]]);inside=center_ids[nearest[,1]] %in% ids
       memberships[[paste0('PT_',j)]]=ifelse(inside,pd$Pseudotime,NA_real_)
       memberships[[paste0('W_',j)]]=as.integer(inside)
       pdefs[[paste0('Path',j)]]=ids
    }
    stopifnot(all(rowSums(memberships[,grep('^W_',names(memberships)),drop=FALSE])>=1))
    fwrite(memberships,file.path(tab,paste0('Monocle2_',name,'_path_cells.csv.gz')))
    write_json(list(root_center=center_root,paths=pdefs,
       assignment='Nearest DDRTree center on each root-to-leaf path; shared trunks appear on multiple paths; no re-fit'),
       file.path(tab,paste0('Monocle2_',name,'_paths.json')),pretty=TRUE,auto_unbox=TRUE)
    saveRDS(res,file.path(output,paste0('Monocle2_',name,'.rds')),compress=FALSE)
    summaries[[name]]=list(requested_root_group=requested,requested_State=state,
       actual_root_cell=root,actual_root_group=as.character(pd[root,'coarse_state']),
       actual_root_State=as.character(pd[root,'State']),root_distance_check=TRUE,
       all_cells=nrow(pd),states=length(unique(pd$State)),geometry_unchanged=TRUE)
 }
 write_json(summaries,file.path(tab,'Monocle2_roots.json'),pretty=TRUE,auto_unbox=TRUE)
}else stop('Unknown method')
write_json(list(method=method,block=key,cells=nrow(m),elapsed_seconds=proc.time()[['elapsed']]-start,
 all_cells_fitted=TRUE,sampled_cells=0,code_file=args[0]),
 file.path(logdir,paste0(key,'_',method,'_completed.json')),pretty=TRUE,auto_unbox=TRUE)
status('complete')
