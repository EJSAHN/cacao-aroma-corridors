"""Prespecified LOCO covariance and PC sensitivity design; no result-based tuning."""
import numpy as np
from ld_audit import family
from adjusted_math import relationship_components,metric_transform

SETTINGS = dict(background_bin_bp=100000,minimum_background_markers=200,
                background_policy='ONE_PER_100KB_BIN_AND_SOURCE_TAG_FAMILY',
                chromosome_policy='LEAVE_TEST_CHROMOSOME_OUT',
                primary_kinship_ridge=1e-6,sensitivity_kinship_ridge=0.01,
                pc_counts=[2,5],kinship_scaling='MEAN_DIAGONAL_ONE',
                intercept_stabilization='PLUS_ONES_MATRIX_OVER_N_REMOVED_BY_GLS_INTERCEPT')
METHODS=('RAW','PC2','PC5','V_RIDGE_1e_minus6','V_RIDGE_001','V_PC2_RIDGE_1e_minus6')

def select_background(markers,held_out_chromosome,bin_bp=100000,minimum=200):
    """Deterministic position-first thinning, performed without any r2 calculation."""
    candidates=[i for i,m in enumerate(markers) if str(m['chromosome'])!=str(held_out_chromosome)]
    candidates.sort(key=lambda i:(int(markers[i]['chromosome']),int(markers[i]['position_bp']),str(markers[i]['marker_key'])))
    selected=[];bins=set();families=set()
    for i in candidates:
        m=markers[i];f=family(m['source_tag_key'])
        if f is None:raise ValueError('Unparseable training source tag family')
        key=(str(m['chromosome']),(int(m['position_bp'])-1)//bin_bp)
        if key in bins or f in families:continue
        selected.append(i);bins.add(key);families.add(f)
    if len(selected)<minimum:raise ValueError('Too few LOCO background markers: '+str(len(selected)))
    return selected

def fit_models(dosages,markers,held_out_chromosome):
    d=np.asarray(dosages,dtype=float)
    if d.ndim!=2 or len(d)!=len(markers) or not np.isfinite(d).all():
        raise ValueError('Complete marker-by-sample dosages are required; no imputation')
    if d.shape[1]<50:raise ValueError('At least 50 samples required')
    selected=select_background(markers,held_out_chromosome,SETTINGS['background_bin_bp'],SETTINGS['minimum_background_markers'])
    r=relationship_components(d[selected],max_pcs=5)
    K=r['K'];pcs=r['pcs'];n=d.shape[1]
    if pcs.shape[1]<5:raise ValueError('Insufficient background rank for PC5 sensitivity')
    intercept_mode=np.ones((n,n))/n
    v=K+SETTINGS['primary_kinship_ridge']*np.eye(n)+intercept_mode
    v2=K+SETTINGS['sensitivity_kinship_ridge']*np.eye(n)+intercept_mode
    specifications=[('RAW',None,None),('PC2',None,pcs[:,:2]),('PC5',None,pcs[:,:5]),
                    ('V_RIDGE_1e_minus6',v,None),('V_RIDGE_001',v2,None),
                    ('V_PC2_RIDGE_1e_minus6',v,pcs[:,:2])]
    models={}
    for name,V,S in specifications:
        T,diagnostics=metric_transform(n,V=V,covariates=S)
        models[name]=dict(T=T,diagnostics=diagnostics)
    return dict(models=models,background_indices=selected,
                background_audit=dict(held_out_chromosome=str(held_out_chromosome),sample_count=n,
                                      available_markers=len(markers),selected_background_markers=len(selected),
                                      **SETTINGS,**r['diagnostics']),
                eigenvalues=r['eigenvalues'])
