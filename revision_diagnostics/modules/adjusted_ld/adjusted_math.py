"""Independent linear-algebra implementation of conditional adjusted dosage r2.

Scientific definitions: Mangin et al. (2012), DOI 10.1038/hdy.2011.73.
For covariance V and nuisance design C=[1,S], the residual precision is
P=V^-1 - V^-1 C (C' V^-1 C)^+ C' V^-1. Squared normalized x' P y
is the generalized partial correlation. The implementation below uses whitening
and an orthogonal nuisance projection, not a translation of LDcorSV R code.
The caller chooses the background markers, covariance regularization and PCs.
Complete calls are required: missing values are rejected, never filled with zero.
"""
from __future__ import annotations
import numpy as np


def _matrix(value, name):
    x=np.asarray(value,dtype=np.float64)
    if x.ndim!=2 or not np.isfinite(x).all():
        raise ValueError(name+' must be a finite two-dimensional matrix')
    return x


def relationship_components(background_dosages,max_pcs=5):
    """Compute a centered dosage relationship kernel and standardized-locus PCs.

    Rows are selected background marker records; columns are samples.
    K is ZZ'/sum(2p(1-p)), normalized to mean diagonal 1. PCs use loci scaled by
    sqrt(2p(1-p)); these are continuous genetic PCs, not ancestry proportions.
    No marker selection, chromosome filtering, shrinkage or imputation occurs here.
    """
    d=_matrix(background_dosages,'background_dosages')
    m,n=d.shape
    if int(max_pcs)!=max_pcs or max_pcs<1:raise ValueError('max_pcs must be a positive integer')
    if m<2 or n<3:raise ValueError('Insufficient background dimensions')
    if not np.isin(d,[0.,1.,2.]).all():raise ValueError('Background must contain 0/1/2 dosages only')
    p=d.mean(axis=1)/2.;heter=2*p*(1-p)
    variable=(heter>0)&(np.var(d,axis=1)>0)
    if int(variable.sum())<max(2,int(max_pcs)):
        raise ValueError('Insufficient variable background markers')
    x=d[variable]-2*p[variable,None]
    denom=float(heter[variable].sum())
    raw=(x.T@x)/denom
    scale=float(np.trace(raw)/n)
    if not np.isfinite(scale) or scale<=0:raise ValueError('Degenerate relationship kernel')
    K=(raw+raw.T)/(2*scale)
    standardized=x/np.sqrt(heter[variable,None])
    G=standardized.T@standardized/len(standardized)
    eig,vec=np.linalg.eigh((G+G.T)/2)
    order=np.argsort(eig)[::-1];eig=eig[order];vec=vec[:,order]
    tol=max(1.,float(eig[0]))*1e-12
    if float(eig[-1]) < -tol:raise ValueError('PCA covariance is not positive semidefinite')
    positive=eig>tol;rank=int(positive.sum())
    if rank<int(max_pcs):raise ValueError('Insufficient PCA rank for requested covariates')
    pcs=vec[:,:int(max_pcs)].copy()
    # Fix signs deterministically for reports; correlations are sign-invariant.
    for j in range(pcs.shape[1]):
        pivot=int(np.argmax(np.abs(pcs[:,j])))
        if pcs[pivot,j]<0:pcs[:,j]*=-1
    evalk=np.linalg.eigvalsh(K)
    kp=evalk[evalk>max(1.,float(evalk[-1]))*1e-12]
    info=dict(background_input_markers=m,background_variable_markers=int(variable.sum()),
              samples=n,vanraden_denominator=denom,raw_mean_diagonal=scale,
              kernel_mean_diagonal=float(np.trace(K)/n),kernel_min_eigenvalue=float(evalk[0]),
              kernel_max_eigenvalue=float(evalk[-1]),kernel_positive_rank=int(len(kp)),
              kernel_positive_condition_number=float(kp[-1]/kp[0]) if len(kp) else None,
              pca_positive_rank=rank,pc_count=int(max_pcs),
              pca_boundary_eigengap=float(eig[int(max_pcs)-1]-eig[int(max_pcs)]) if len(eig)>int(max_pcs) else None,
              pca_variance_fraction=[float(v/eig[positive].sum()) for v in eig[:int(max_pcs)]],
              missing_values_imputed=False)
    return dict(K=K,pcs=pcs,diagnostics=info,eigenvalues=eig)


def metric_transform(n_samples,V=None,covariates=None,eigen_rtol=1e-12):
    """Return T with P=T T', where P is nuisance-adjusted inverse covariance.

    With V=None this is ordinary/PC residual correlation. With V, use its positive
    eigenspace (Moore-Penrose metric), then remove the intercept and any supplied
    covariates in whitened space. A covariance whose retained precision annihilates
    the intercept is rejected; callers must explicitly regularize centered K.
    """
    n=int(n_samples)
    if n<3 or n!=n_samples:raise ValueError('At least three samples required')
    if eigen_rtol<=0 or eigen_rtol>=1:raise ValueError('Invalid eigen tolerance')
    if V is None:
        white=np.eye(n);values=np.ones(n)
    else:
        v=_matrix(V,'V')
        if v.shape!=(n,n) or not np.allclose(v,v.T,rtol=1e-12,atol=1e-12):
            raise ValueError('V must be symmetric and sample-aligned')
        values,vectors=np.linalg.eigh((v+v.T)/2)
        cutoff=max(float(np.max(np.abs(values))),1.)*eigen_rtol
        if values[0]<-cutoff:raise ValueError('V is not positive semidefinite')
        keep=values>cutoff
        if int(keep.sum())<3:raise ValueError('Covariance retained rank is too small')
        values=values[keep]
        white=vectors[:,keep]/np.sqrt(values)[None,:]
    cov=np.empty((n,0)) if covariates is None else _matrix(covariates,'covariates')
    if cov.shape[0]!=n:raise ValueError('Covariates are not sample-aligned')
    design=np.column_stack((np.ones(n),cov))
    weighted=white.T@design
    if np.linalg.norm(weighted[:,0])<=1e-10:
        raise ValueError('Intercept is not estimable in covariance metric; explicitly regularize V')
    u,s,_=np.linalg.svd(weighted,full_matrices=True)
    cutoff=max(float(s[0]),1.)*eigen_rtol
    rank=int(np.sum(s>cutoff))
    if rank!=design.shape[1]:raise ValueError('Nuisance design is rank deficient')
    residual_dimension=white.shape[1]-rank
    if residual_dimension<2:raise ValueError('Insufficient residual dimension')
    transform=white@u[:,rank:]
    ortho=float(np.max(np.abs(design.T@transform)))
    if ortho>1e-7*max(1.,float(np.linalg.norm(transform))):
        raise ValueError('Nuisance projection failed orthogonality check')
    info=dict(samples=n,covariance_retained_rank=int(white.shape[1]),
              covariance_min_retained_eigenvalue=float(values.min()),
              covariance_max_eigenvalue=float(values.max()),
              covariance_retained_condition_number=float(values.max()/values.min()),
              nuisance_rank=rank,residual_dimension=residual_dimension,
              nuisance_orthogonality_max_abs=ortho,eigen_relative_tolerance=eigen_rtol)
    return transform,info


def project_dosages(dosages,transform):
    """Project marker rows and retain scale-aware zero-residual-variance flags."""
    d=_matrix(dosages,'dosages');t=_matrix(transform,'transform')
    if d.shape[1]!=t.shape[0]:raise ValueError('Dosage/transform sample dimensions differ')
    projected=d@t
    ss=np.einsum('ij,ij->i',projected,projected)
    centered=d-d.mean(axis=1,keepdims=True)
    rawss=np.einsum('ij,ij->i',centered,centered)
    valid=ss>1e-12*np.maximum(1.,rawss)
    unit=np.zeros_like(projected)
    unit[valid]=projected[valid]/np.sqrt(ss[valid,None])
    return dict(unit=unit,valid=valid,ss=ss,samples=d.shape[1])


def correlations_from_projected(a,b,min_samples=50):
    """Compute matching-row correlations after an identical metric transform."""
    if a['unit'].shape!=b['unit'].shape or a['samples']!=b['samples']:
        raise ValueError('Projected pairs are not sample/dimension aligned')
    size=len(a['unit']);n=int(a['samples'])
    enough=n>=int(min_samples)
    valid=a['valid']&b['valid']&enough
    r2=np.full(size,np.nan)
    dot=np.einsum('ij,ij->i',a['unit'][valid],b['unit'][valid])
    if np.any(np.abs(dot)>1+1e-10):raise ValueError('Adjusted correlation is out of range')
    r2[valid]=np.clip(dot*dot,0,1)
    return dict(r2=r2,n=np.full(size,n,dtype=np.int32),
                status=np.where(not enough,1,np.where(valid,0,2)).astype(np.uint8))


def pair_correlations(xrows,yrows,transform,min_samples=50):
    """Matching rows, same sample order: returns r2/n/status vectors.

    Status 0=valid, 1=insufficient samples, 2=zero residual variance. Missing
    source values cause an error before projection; they are never imputed.
    """
    x=_matrix(xrows,'xrows');y=_matrix(yrows,'yrows')
    if x.shape!=y.shape:raise ValueError('Pair rows must have identical dimensions')
    return correlations_from_projected(project_dosages(x,transform),project_dosages(y,transform),min_samples)


def metric_pairvalues(dosages,pairs,transform,min_samples=50,batch=4096):
    """Project each locus once, then score indexed pairs in bounded batches."""
    projected=project_dosages(dosages,transform)
    pp=np.asarray(pairs,dtype=np.int64)
    if pp.ndim!=2 or pp.shape[1]!=2:raise ValueError('Pairs must have two columns')
    if len(pp) and (pp.min()<0 or pp.max()>=len(projected['unit']) or np.any(pp[:,0]>=pp[:,1])):
        raise ValueError('Invalid unordered marker pair')
    if len(pp) and len(np.unique(pp[:,0]*len(projected['unit'])+pp[:,1]))!=len(pp):
        raise ValueError('Duplicate marker pairs')
    if batch<1:raise ValueError('Invalid batch size')
    out=dict(r2=np.full(len(pp),np.nan),n=np.zeros(len(pp),np.int32),status=np.ones(len(pp),np.uint8))
    for start in range(0,len(pp),batch):
        stop=min(start+batch,len(pp));a,b=pp[start:stop].T
        pa={k:v[a] if isinstance(v,np.ndarray) else v for k,v in projected.items()}
        pb={k:v[b] if isinstance(v,np.ndarray) else v for k,v in projected.items()}
        value=correlations_from_projected(pa,pb,min_samples)
        for k in out:out[k][start:stop]=value[k]
    return out
