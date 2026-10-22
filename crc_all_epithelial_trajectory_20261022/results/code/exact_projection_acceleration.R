# Source-only definition. Installation affects only this R process, not packages.
make_exact_projection_accelerator <- function(root, min_vertices=1000L, audit_file=NULL) {
  private <- new.env(parent=baseenv())
  private$reference <- princurve::project_to_curve
  private$audit <- list()
  cache=file.path(root,'.runtime/exact_projection_cache')
  dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  Rcpp::sourceCpp(file.path(root,'results/code/exact_projection_bvh.cpp'),
                  env=private,cacheDir=cache,showOutput=FALSE)
  wrapper <- function(x,s,stretch=2) {
    if(is.matrix(s) && nrow(s)>=min_vertices && isTRUE(stretch==0) &&
       ncol(x)==ncol(s) && ncol(x)>0 && nrow(x)>0 &&
       all(is.finite(s)) && all(is.finite(x))) {
      if(!is.null(audit_file)){
        cat(format(Sys.time()),'exact segment search start:',nrow(x),'rows,',nrow(s),'vertices\n')
        flush.console()
      }
      result=private$project_all_rows_bvh(x,s)
      if(isTRUE(result$usable)) {
        audit=result[setdiff(names(result),c('projection','usable'))]
        private$audit[[length(private$audit)+1L]]=audit
        if(!is.null(audit_file)){
          jsonlite::write_json(private$audit,audit_file,pretty=TRUE,auto_unbox=TRUE,digits=NA)
          cat(format(Sys.time()),'exact segment search complete:',
              audit$evaluated_segment_pairs,'of',audit$full_enumeration_pairs,'pairs evaluated\n')
          flush.console()
        }
        return(result$projection)
      }
    }
    private$reference(x=x,s=s,stretch=stretch)
  }
  locate=function(){
    current=environment(methods::selectMethod('getCurves','PseudotimeOrdering'))
    while(!identical(current,emptyenv()) &&
          !exists('project_to_curve',envir=current,inherits=FALSE))current=parent.env(current)
    stopifnot(!identical(current,emptyenv()));current
  }
  replace=function(value){
    current=locate();locked=bindingIsLocked('project_to_curve',current)
    if(locked)unlockBinding('project_to_curve',current)
    assign('project_to_curve',value,envir=current)
    if(locked)lockBinding('project_to_curve',current)
    invisible(environmentName(current))
  }
  list(project=wrapper,reference=private$reference,
    audit=function()private$audit,
    install=function(){
      stopifnot(identical(get('project_to_curve',envir=locate(),inherits=FALSE),private$reference))
      replace(wrapper)
    },uninstall=function()replace(private$reference))
}
