# Candidate runtime-only acceleration. Native princurve projection remains the
# reference kernel; all X rows are always projected. No installed files change.
# This source alone does not install the wrapper or start a biological fit.
make_polyline_projection_accelerator <- function(min_vertices=1000L, tolerance=1e-12,
                                                 max_vertices=256L) {
  private <- new.env(parent=baseenv())
  private$reference <- princurve::project_to_curve
  private$audit <- list()
  Rcpp::cppFunction(env=private, code='
Rcpp::List prune_verified_polyline(Rcpp::NumericMatrix s, double tolerance,
                                   int max_vertices) {
  const int nr=s.nrow(), nc=s.ncol();
  std::vector<int> unique, keep;
  unique.push_back(0);
  for(int i=0;i<nr;++i) for(int k=0;k<nc;++k)
    if(!R_finite(s(i,k))) return Rcpp::List::create(Rcpp::Named("safe")=false);
  for(int i=1;i<nr;++i){
    bool same=true;
    for(int k=0;k<nc;++k) if(s(i,k)!=s(unique.back(),k)){same=false;break;}
    if(!same) unique.push_back(i);
  }
  if(unique.size()<2) return Rcpp::List::create(Rcpp::Named("safe")=false);
  unique.back()=nr-1;
  keep.push_back(unique.front());
  const double tol2=tolerance*tolerance;
  for(size_t j=2;j<unique.size();++j){
    int a=keep.back(), b=unique[j-1], c=unique[j];
    double denominator=0, numerator=0;
    for(int k=0;k<nc;++k){
      double ac=s(c,k)-s(a,k), ab=s(b,k)-s(a,k);
      denominator+=ac*ac; numerator+=ab*ac;
    }
    bool collinear=false;
    if(denominator>0 && R_finite(denominator)){
      double t=numerator/denominator, error2=0;
      for(int k=0;k<nc;++k){
        double error=s(b,k)-s(a,k)-t*(s(c,k)-s(a,k));
        error2+=error*error;
      }
      collinear=t>=0 && t<=1 && error2<=tol2;
    }
    if(!collinear) keep.push_back(b);
    if(keep.size()>static_cast<size_t>(max_vertices))
      return Rcpp::List::create(Rcpp::Named("safe")=false);
  }
  keep.push_back(unique.back());
  if(keep.size()>static_cast<size_t>(max_vertices) || keep.size()>=static_cast<size_t>(nr))
    return Rcpp::List::create(Rcpp::Named("safe")=false);
  // Verify EVERY original vertex against its enclosing retained chord, including
  // monotone position along that chord. This rejects accumulated curvature.
  double max_error2=0;
  for(size_t j=1;j<keep.size();++j){
    int a=keep[j-1], c=keep[j]; double denominator=0;
    for(int k=0;k<nc;++k){double v=s(c,k)-s(a,k);denominator+=v*v;}
    if(!(denominator>0) || !R_finite(denominator))
      return Rcpp::List::create(Rcpp::Named("safe")=false);
    double previous_t=0;
    for(int b=a+1;b<c;++b){
      double numerator=0,error2=0;
      for(int k=0;k<nc;++k) numerator+=(s(b,k)-s(a,k))*(s(c,k)-s(a,k));
      double t=numerator/denominator;
      if(t < -1e-14 || t>1+1e-14 || t+1e-14<previous_t)
        return Rcpp::List::create(Rcpp::Named("safe")=false);
      previous_t=t;
      for(int k=0;k<nc;++k){
        double error=s(b,k)-s(a,k)-t*(s(c,k)-s(a,k));error2+=error*error;
      }
      if(!R_finite(error2) || error2>tol2)
        return Rcpp::List::create(Rcpp::Named("safe")=false);
      if(error2>max_error2)max_error2=error2;
    }
  }
  Rcpp::IntegerVector index(keep.size());
  for(size_t j=0;j<keep.size();++j)index[j]=keep[j]+1;
  return Rcpp::List::create(Rcpp::Named("safe")=true,
    Rcpp::Named("index")=index,Rcpp::Named("max_vertex_error")=std::sqrt(max_error2));
}')
  wrapper <- function(x, s, stretch=2) {
    if(is.matrix(s) && nrow(s)>=min_vertices && isTRUE(stretch==0)) {
      checked <- private$prune_verified_polyline(s, tolerance, max_vertices)
      if(isTRUE(checked$safe)) {
        original_vertices <- nrow(s)
        s <- s[checked$index,,drop=FALSE]
        private$audit[[length(private$audit)+1L]] <- list(
          projected_rows=nrow(x), dimensions=ncol(x), original_vertices=original_vertices,
          retained_vertices=nrow(s), max_vertex_error=checked$max_vertex_error,
          tolerance=tolerance, all_original_rows_projected=TRUE)
      }
    }
    private$reference(x=x,s=s,stretch=stretch)
  }
  list(project=wrapper, reference=private$reference,
       audit=function()private$audit, prune=private$prune_verified_polyline,
       install=function(){
         current=environment(methods::selectMethod('getCurves','PseudotimeOrdering'))
         while(!identical(current,emptyenv()) &&
               !exists('project_to_curve',envir=current,inherits=FALSE))current=parent.env(current)
         stopifnot(!identical(current,emptyenv()),
                   identical(get('project_to_curve',envir=current,inherits=FALSE),private$reference))
         locked=bindingIsLocked('project_to_curve',current)
         if(locked)unlockBinding('project_to_curve',current)
         assign('project_to_curve',wrapper,envir=current)
         if(locked)lockBinding('project_to_curve',current)
         invisible(environmentName(current))
       },
       uninstall=function(){
         current=environment(methods::selectMethod('getCurves','PseudotimeOrdering'))
         while(!identical(current,emptyenv()) &&
               !exists('project_to_curve',envir=current,inherits=FALSE))current=parent.env(current)
         stopifnot(!identical(current,emptyenv()))
         locked=bindingIsLocked('project_to_curve',current)
         if(locked)unlockBinding('project_to_curve',current)
         assign('project_to_curve',private$reference,envir=current)
         if(locked)lockBinding('project_to_curve',current)
       })
}
