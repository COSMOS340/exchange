"""All-gene clustering metrics; no cell or gene subsampling for silhouette."""
import numpy as np
from scipy.cluster.vq import ClusterError, kmeans2
from scipy.spatial.distance import cdist


def standardize_profiles(curves):
    curves=np.asarray(curves,dtype=float)
    mean=curves.mean(axis=1,keepdims=True)
    sd=curves.std(axis=1,keepdims=True)
    return np.divide(curves-mean,sd,out=np.zeros_like(curves),where=sd>1e-12)


def adjusted_rand(a,b):
    a=np.asarray(a);b=np.asarray(b)
    if len(a)!=len(b):raise ValueError('Partitions differ in length')
    if len(a)<2:return 1.
    _,a=np.unique(a,return_inverse=True);_,b=np.unique(b,return_inverse=True)
    contingency=np.zeros((a.max()+1,b.max()+1),dtype=np.int64)
    np.add.at(contingency,(a,b),1)
    choose=lambda x: np.asarray(x,dtype=float)*(np.asarray(x,dtype=float)-1)/2
    observed=choose(contingency).sum()
    sa,sb=choose(contingency.sum(axis=1)).sum(),choose(contingency.sum(axis=0)).sum()
    expected=sa*sb/choose(len(a))
    denominator=(sa+sb)/2-expected
    return 1. if abs(denominator)<1e-15 else float((observed-expected)/denominator)


def cluster_all_profiles(x,k,seed):
    x=np.asarray(x,dtype=float)
    if not np.isfinite(x).all() or len(x)<=k:
        raise ValueError('Need more finite gene profiles than modules')
    best=None
    for start in range(5):
        try:
            centers,labels=kmeans2(x,k,iter=100,minit='++',missing='raise',
                                  rng=np.random.default_rng(seed+start))
            if len(np.unique(labels))!=k:continue
            inertia=float(np.square(x-centers[labels]).sum())
            if best is None or inertia<best['inertia']:
                best=dict(centers=centers,labels=labels,inertia=inertia,initialization=start)
        except (ValueError,ClusterError):
            continue
    if best is None:raise ValueError('Every initialization produced an empty module')
    return best


def exact_silhouettes(x,partitions,chunk=128):
    """One blocked distance pass evaluates every candidate K on every gene."""
    keys=sorted(partitions)
    sizes={k:np.bincount(partitions[k],minlength=k) for k in keys}
    indicators=np.concatenate([np.eye(k)[partitions[k]] for k in keys],axis=1)
    scores={k:np.zeros(len(x)) for k in keys}
    for lo in range(0,len(x),chunk):
        hi=min(lo+chunk,len(x))
        sums=cdist(x[lo:hi],x)@indicators
        offset=0
        for k in keys:
            group=partitions[k][lo:hi]
            counts=sizes[k][group]
            means=sums[:,offset:offset+k]/sizes[k][None,:]
            within_sum=sums[np.arange(hi-lo),offset+group]
            within=np.divide(within_sum,counts-1,out=np.zeros(hi-lo),where=counts>1)
            means[np.arange(hi-lo),group]=np.inf
            between=means.min(axis=1)
            value=np.divide(between-within,np.maximum(within,between),
                            out=np.zeros_like(within),where=np.maximum(within,between)>0)
            value[counts==1]=0
            scores[k][lo:hi]=value
            offset+=k
    return scores
