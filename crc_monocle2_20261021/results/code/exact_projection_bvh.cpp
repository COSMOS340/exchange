// GPL-3. Projection arithmetic/output ordering adapted from princurve 2.1.6
// (https://github.com/rcannood/princurve, src/project_to_curve.cpp).
// All original curve segments and all input rows are retained. The hierarchy
// only skips boxes whose conservative distance cannot beat the current answer.
#include <Rcpp.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <functional>
#include <limits>
#include <utility>
#include <vector>
using namespace Rcpp;

struct ProjectionBox {
  int begin, end, left=-1, right=-1;
  std::vector<double> lower, upper;
  ProjectionBox(int a,int b,int p):begin(a),end(b),lower(p),upper(p){}
};

// [[Rcpp::export]]
List project_all_rows_bvh(NumericMatrix x, NumericMatrix s) {
  const int ns=s.nrow()-1, nx=x.nrow(), p=x.ncol();
  if(ns<1 || nx<1 || p<1 || s.ncol()!=p)stop("Invalid projection dimensions");
  double scale=1;
  for(int k=0;k<p;++k){
    for(int j=0;j<=ns;++j){
      if(!R_finite(s(j,k)))stop("Nonfinite curve coordinate");
      scale=std::max(scale,std::abs(s(j,k)));
    }
    for(int i=0;i<nx;++i){
      if(!R_finite(x(i,k)))stop("Nonfinite point coordinate");
      scale=std::max(scale,std::abs(x(i,k)));
    }
  }
  NumericMatrix delta=no_init(ns,p);
  NumericVector lengths=no_init(ns);
  bool any_nondegenerate=false;
  for(int j=0;j<ns;++j){
    double value=0;
    for(int k=0;k<p;++k){
      double d=s(j+1,k)-s(j,k);delta[k*ns+j]=d;value+=d*d;
    }
    lengths[j]=value;
    any_nondegenerate|=(value>0 && R_finite(value));
  }
  if(!any_nondegenerate)return List::create(_["usable"]=false);
  std::vector<ProjectionBox> nodes;
  nodes.reserve(2*ns);
  std::function<int(int,int)> build=[&](int first,int last){
    int id=nodes.size();nodes.emplace_back(first,last,p);
    if(last-first<=16){
      for(int k=0;k<p;++k){
        double low=s(first,k),high=low;
        for(int j=first+1;j<=last;++j){
          low=std::min(low,s(j,k));high=std::max(high,s(j,k));
        }
        nodes[id].lower[k]=low;nodes[id].upper[k]=high;
      }
    }else{
      int middle=first+(last-first)/2;
      int left=build(first,middle),right=build(middle,last);
      nodes[id].left=left;nodes[id].right=right;
      for(int k=0;k<p;++k){
        nodes[id].lower[k]=std::min(nodes[left].lower[k],nodes[right].lower[k]);
        nodes[id].upper[k]=std::max(nodes[left].upper[k],nodes[right].upper[k]);
      }
    }
    return id;
  };
  int root=build(0,ns);
  NumericMatrix output=no_init(nx,p);
  NumericVector lambda=no_init(nx),distance=no_init(nx);
  NumericVector point=no_init(p),candidate=no_init(p),best_point=no_init(p);
  // The generous scale-dependent guard covers floating coordinate projection,
  // subtraction, squares and p-term summation. Boxes use long-double bounds.
  const long double guard=4096.L*DBL_EPSILON*p*p*scale*scale;
  auto lower_bound=[&](int id){
    long double value=0;
    for(int k=0;k<p;++k){
      long double gap=0;
      if(point[k]<nodes[id].lower[k])gap=(long double)nodes[id].lower[k]-point[k];
      if(point[k]>nodes[id].upper[k])gap=(long double)point[k]-nodes[id].upper[k];
      value+=gap*gap;
    }
    return value;
  };
  std::vector<std::pair<int,long double>> stack;stack.reserve(128);
  double evaluated=0;
  for(int i=0;i<nx;++i){
    if((i&1023)==0)checkUserInterrupt();
    for(int k=0;k<p;++k)point[k]=x(i,k);
    double best_distance=R_PosInf,best_lambda=-1;
    int best_segment=std::numeric_limits<int>::max();
    stack.clear();stack.emplace_back(root,0);
    while(!stack.empty()){
      auto item=stack.back();stack.pop_back();
      if(item.second>(long double)best_distance+guard)continue;
      const ProjectionBox& node=nodes[item.first];
      if(node.left>=0){
        long double left=lower_bound(node.left),right=lower_bound(node.right);
        if(left<=right){stack.emplace_back(node.right,right);stack.emplace_back(node.left,left);}
        else{stack.emplace_back(node.left,left);stack.emplace_back(node.right,right);}
        continue;
      }
      for(int j=node.begin;j<node.end;++j){
        ++evaluated;
        // Same scalar operations, order and zero-length behavior as reference.
        double v=0;
        for(int k=0;k<p;++k)v+=delta(j,k)*(point[k]-s(j,k));
        v/=lengths(j);
        if(v<0)v=0.0;if(v>1)v=1.0;
        double w=0;
        for(int k=0;k<p;++k){
          double u=s(j,k)+v*delta(j,k);candidate[k]=u;
          w+=(u-point[k])*(u-point[k]);
        }
        // Native enumeration chooses the earliest segment on an exact tie.
        if(w<best_distance || (w==best_distance && R_finite(w) && j<best_segment)){
          best_distance=w;best_segment=j;best_lambda=j+.1+.9*v;
          for(int k=0;k<p;++k)best_point[k]=candidate[k];
        }
      }
    }
    if(!R_finite(best_distance))return List::create(_["usable"]=false);
    lambda[i]=best_lambda;distance[i]=best_distance;
    for(int k=0;k<p;++k)output[k*nx+i]=best_point[k];
  }
  // Match the native pair<int,double> sort and subsequent chord-length sum.
  std::vector<std::pair<int,double>> keys;keys.reserve(nx);
  for(int i=0;i<nx;++i)keys.emplace_back(i,lambda[i]);
  std::sort(keys.begin(),keys.end(),[](const std::pair<int,double>& a,
                                    const std::pair<int,double>& b){return a.second<b.second;});
  IntegerVector ordering=no_init(nx);
  for(int i=0;i<nx;++i)ordering[i]=keys[i].first;
  lambda[ordering[0]]=0;
  for(int i=1;i<nx;++i){
    int a=ordering[i-1],b=ordering[i];double w=0;
    for(int k=0;k<p;++k){double v=output(b,k)-output(a,k);w+=v*v;}
    lambda[b]=lambda[a]+std::sqrt(w);
  }
  output.attr("dimnames")=List::create(rownames(x),colnames(x));
  distance.attr("names")=rownames(x);lambda.attr("names")=rownames(x);
  List projection;
  projection["s"]=output;projection["ord"]=ordering+1;
  projection["lambda"]=lambda;projection["dist_ind"]=distance;
  double total_distance=sum(distance);
  projection["dist"]=total_distance;projection.attr("class")="principal_curve";
  return List::create(_["usable"]=true,_["projection"]=projection,
    _["evaluated_segment_pairs"]=evaluated,
    _["full_enumeration_pairs"]=(double)nx*ns,
    _["curve_vertices_retained"]=s.nrow(),_["input_rows_retained"]=nx,
    _["rounding_guard"]=(double)guard);
}
