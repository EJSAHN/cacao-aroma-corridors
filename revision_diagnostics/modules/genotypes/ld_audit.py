"""Conditional dosage correlations, source-tag sensitivity and peak coordinate anchors.
A/M/C = 0/1/2 is an explicit unproven interpretation. No structure adjustment.
"""
from collections import Counter, defaultdict
import re
import numpy as np
SCOPE='CONDITIONAL_A_M_C_DOSAGE;UNADJUSTED;OBSERVED_COORDINATE_ANCHORS_ONLY'

def scalar(v):
    if isinstance(v,(np.integer,)): return int(v)
    if isinstance(v,(np.bool_,)): return bool(v)
    if isinstance(v,(np.floating,float)): return float(v) if np.isfinite(v) else None
    return v

def family(key):
    m=re.fullmatch(r'(.*)--([0-9]+)',str(key or ''))
    return m.group(1) if m else None

def local_pairs(markers,maximum,limit=1500000):
    by=defaultdict(list)
    for i,m in enumerate(markers): by[str(m['chromosome'])].append(i)
    chunks=[];count=0
    for chrom in sorted(by):
        ids=np.array(sorted(by[chrom],key=lambda i:(int(markers[i]['position_bp']),i)),dtype=np.int64)
        pos=np.array([int(markers[i]['position_bp']) for i in ids]);ends=np.searchsorted(pos,pos+maximum,side='right')
        count+=int(np.sum(ends-np.arange(len(ids))-1))
        if count>limit: raise ValueError('Local pair census exceeds safety limit; not sampled')
        for k,end in enumerate(ends):
            others=ids[k+1:int(end)]
            if len(others):chunks.append(np.column_stack((np.minimum(ids[k],others),np.maximum(ids[k],others))))
    if not chunks:return np.empty((0,2),dtype=np.int64)
    p=np.concatenate(chunks);return p[np.lexsort((p[:,1],p[:,0]))]

def correlation(dosages,pairs,minimum=50,batch=2048):
    d,p=np.asarray(dosages),np.asarray(pairs,dtype=np.int64)
    if d.ndim!=2 or p.ndim!=2 or p.shape[1]!=2:raise ValueError('Invalid LD dimensions')
    if np.any(np.isinf(d)) or not np.all(np.isin(d[np.isfinite(d)],[0,1,2])):raise ValueError('Invalid dosage')
    if len(p) and (p.min()<0 or p.max()>=len(d) or np.any(p[:,0]>=p[:,1])):raise ValueError('Invalid LD pair')
    if len(p) and len(np.unique(p[:,0]*len(d)+p[:,1]))!=len(p):raise ValueError('Duplicate LD pairs')
    n=np.zeros(len(p),dtype=np.int32);r2=np.full(len(p),np.nan);status=np.ones(len(p),dtype=np.uint8)
    for lo in range(0,len(p),batch):
        hi=min(lo+batch,len(p));q=p[lo:hi];x=d[q[:,0]].astype(float);y=d[q[:,1]].astype(float)
        present=np.isfinite(x)&np.isfinite(y);x=np.where(present,x,0);y=np.where(present,y,0)
        nn=present.sum(axis=1);sx=x.sum(axis=1);sy=y.sum(axis=1)
        cov=nn*(x*y).sum(axis=1)-sx*sy;vx=nn*(x*x).sum(axis=1)-sx*sx;vy=nn*(y*y).sum(axis=1)-sy*sy
        enough=nn>=minimum;valid=enough&(vx>0)&(vy>0);rr=np.full(hi-lo,np.nan)
        rr[valid]=cov[valid]**2/(vx[valid]*vy[valid])
        if np.any(rr[valid]>1+1e-12):raise ValueError('Correlation out of bounds')
        n[lo:hi]=nn;r2[lo:hi]=np.clip(rr,0,1);status[lo:hi]=np.where(~enough,1,np.where(valid,0,2))
    return dict(n=n,r2=r2,status=status)

def statistics(values,mask):
    rr,nn,st=(values[k][mask] for k in ('r2','n','status'));v=rr[np.isfinite(rr)]
    r=dict(pair_count=len(rr),valid_pairs=len(v),insufficient_samples=int(np.sum(st==1)),zero_pair_variance=int(np.sum(st==2)),
           mean_r2=float(v.mean()) if len(v) else None,median_r2=float(np.median(v)) if len(v) else None,
           q90_r2=float(np.quantile(v,.9)) if len(v) else None,q95_r2=float(np.quantile(v,.95)) if len(v) else None,
           min_complete_samples=int(nn.min()) if len(nn) else None,max_complete_samples=int(nn.max()) if len(nn) else None)
    for t in (.2,.5,.8):
        k=str(t).replace('.','_');r['pairs_r2_ge_'+k]=int(np.sum(v>=t));r['fraction_r2_ge_'+k]=float(np.mean(v>=t)) if len(v) else None
    return r

def equal(a,b):
    a,b=scalar(a),scalar(b)
    return (a is None and b is None) or (a is not None and b is not None and
           (abs(float(a)-float(b))<=1e-12 if isinstance(a,(int,float)) and isinstance(b,(int,float)) else a==b))

def compare(ctx,label,observed,expected,fields):
    diff=[dict(comparison=label,field=k,observed=scalar(observed.get(k)),frozen=scalar(expected.get(k)))
          for k in fields if not equal(observed.get(k),expected.get(k))]
    ctx.check(label,len(diff),0);return diff

def neighbors(panel,pairs,dist,values,include,window=50000):
    n=len(panel.markers);local=np.flatnonzero((dist>0)&(dist<=window)&include)
    eligible=np.bincount(pairs[local].ravel(),minlength=n);valid=local[np.isfinite(values['r2'][local])]
    evaluable=np.bincount(pairs[valid].ravel(),minlength=n);best=np.full(n,np.nan);partner=np.full(n,-1,dtype=int);distance=np.full(n,-1,dtype=int)
    for z in valid:
        i,j=map(int,pairs[z]);r=float(values['r2'][z]);d=int(dist[z])
        for k,other in ((i,j),(j,i)):
            if partner[k]<0 or r>best[k] or (r==best[k] and other<partner[k]):best[k],partner[k],distance[k]=r,other,d
    rows=[]
    for i,m in enumerate(panel.markers):
        j=int(partner[i]);rows.append(dict(marker_key=m['marker_key'],chromosome=str(m['chromosome']),position_bp=int(m['position_bp']),window_bp=window,
          eligible_other_markers=int(eligible[i]),valid_other_markers=int(evaluable[i]),best_partner_key=panel.markers[j]['marker_key'] if j>=0 else '',
          best_partner_distance_bp=int(distance[i]) if j>=0 else None,best_r2=scalar(best[i]),
          status='EVALUABLE' if j>=0 else 'NO_VALID_PAIR' if eligible[i] else 'NO_OTHER_MARKER_IN_WINDOW'))
    summary=dict(window_bp=window,typed_markers=n,markers_with_other_marker=int(np.sum(eligible>0)),markers_with_valid_pair=int(np.sum(evaluable>0)),
                 markers_without_any_neighbor=int(np.sum(eligible==0)),markers_with_neighbors_but_no_valid_pair=int(np.sum((eligible>0)&(evaluable==0))))
    for t in (.2,.5,.8):
        key=str(t).replace('.','_');nt=int(np.sum(best>=t));ne=int(np.sum(evaluable>0))
        summary['markers_best_r2_ge_'+key]=nt;summary['fraction_all_markers_ge_'+key]=nt/n if n else None;summary['fraction_evaluable_ge_'+key]=nt/ne if ne else None
    return rows,summary

def scalar_self_checks(ctx):
    rng=np.random.default_rng(20261026);d=rng.integers(0,3,(16,73)).astype(float);d[rng.random(d.shape)<.1]=np.nan
    pairs=np.array([(i,j) for i in range(len(d)) for j in range(i+1,len(d))]);v=correlation(d,pairs,10,7);errors=[]
    for z,(i,j) in enumerate(pairs):
        ok=np.isfinite(d[i])&np.isfinite(d[j]);errors.append(abs(v['r2'][z]-np.corrcoef(d[i,ok],d[j,ok])[0,1]**2))
    ctx.check('LD independent scalar Pearson reference',bool(max(errors)<1e-12),True)
    s=correlation(np.array([[0,1,2,np.nan],[0,1,2,1],[1,1,1,1.]]),np.array([[0,1],[0,2],[1,2]]),3)
    ctx.check('LD missing intersection and zero variance',(s['n'].tolist(),s['status'].tolist()),([3,3,4],[0,2,2]))
    ctx.check('LD perfect correlation',scalar(s['r2'][0]),1.0);ctx.check('LD missing is not zero',bool(np.isnan(s['r2'][1:]).all()),True)

def peak_tables(ctx,panels,lookup):
    peaks=ctx.read('source_data/02_association_units.xlsx','unique_peaks')
    ledger=ctx.read('alignment_coverage/02_fixed_marker_cohorts.xlsx','marker_record_ledger');positions=defaultdict(list);exact_records=Counter()
    for r in ledger:
        if r['in_C'] and r['resource'] in ('Diversity','Nacional_G7','Nacional_MAF5'):
            key=(r['resource'],str(r['reference_chromosome']));bp=int(r['reference_snp_bp']);positions[key].append(bp);exact_records[key+(bp,)]+=1
    positions={k:np.array(sorted(set(v)),dtype=np.int64) for k,v in positions.items()}
    anchors={name:{(str(m['chromosome']),int(m['position_bp'])):i for i,m in enumerate(p.markers)} for name,p in panels.items()}
    output=[];counts=Counter()
    for peak in peaks:
        chrom,bp=str(peak['chromosome']),int(peak['start_bp'])
        if int(peak['end_bp'])!=bp:raise ValueError('Expected point peak')
        for resource in ('Diversity','Nacional_G7','Nacional_MAF5'):
            pp=positions.get((resource,chrom),np.empty(0,dtype=np.int64));p=int(np.searchsorted(pp,bp))
            nearest=min([abs(int(pp[i])-bp) for i in (p-1,p) if 0<=i<len(pp)],default=None)
            exact=exact_records[(resource,chrom,bp)];number50=int(np.searchsorted(pp,bp+50000,side='right')-np.searchsorted(pp,bp-50000,side='left'))
            anchor=anchors.get(resource,{}).get((chrom,bp));groups=panels[resource].groups if resource in panels else {'NOT_APPLICABLE':[]}
            for group in groups:
                policies=('ALL_LOCAL_PAIRS','EXCLUDE_SAME_SOURCE_TAG_FAMILY') if resource in panels else ('NOT_APPLICABLE',)
                for policy in policies:
                    nr=lookup.get((resource,group,policy,anchor)) if anchor is not None else None
                    if resource=='Diversity':status='EXACT_COORDINATE_MATCH_NO_ACCESSION_GENOTYPES' if exact else 'NEARBY_ONLY_NO_ACCESSION_GENOTYPES' if number50 else 'NO_MARKER_WITHIN_50KB_NO_ACCESSION_GENOTYPES'
                    elif anchor is not None:
                        if nr['status']!='EVALUABLE': status='EXACT_COORDINATE_ANCHOR_NO_EVALUABLE_NEIGHBOR'
                        elif nr['best_r2']>=.8: status='EXACT_COORDINATE_ANCHOR_WITH_OBSERVED_HIGH_CORRELATION_NEIGHBOR'
                        else: status='EXACT_COORDINATE_ANCHOR_EVALUABLE_NEIGHBORS_BELOW_R2_0_8'
                    elif exact:status='EXACT_COORDINATE_MATCH_EXCLUDED_BY_FIXED_LD_QC'
                    elif number50:status='NEARBY_ONLY_FOCAL_GENOTYPE_UNAVAILABLE'
                    else:status='NO_MARKER_WITHIN_50KB_FOCAL_GENOTYPE_UNAVAILABLE'
                    row=dict(query_id=peak['query_id'],chromosome=chrom,position_bp=bp,n_events=peak['n_events'],traits=peak['traits'],resource=resource,sample_group=group,pair_policy=policy,
                      corrected_exact_records=exact,corrected_exact_coordinate_match=bool(exact),nearest_corrected_marker_distance_bp=nearest,
                      corrected_unique_positions_within_50kb=number50,has_ld_qc_coordinate_anchor=anchor is not None,
                      coordinate_anchor_marker_key=panels[resource].markers[anchor]['marker_key'] if anchor is not None else '',
                      ld_assessment_status=nr['status'] if nr else 'NOT_ASSESSABLE_FOCAL_GENOTYPE_UNAVAILABLE' if resource!='Diversity' else 'NOT_ASSESSABLE_NO_ACCESSION_GENOTYPES',
                      valid_other_qc_markers_within_50kb=nr['valid_other_markers'] if nr else None,best_observed_partner_key=nr['best_partner_key'] if nr else '',
                      best_observed_partner_distance_bp=nr['best_partner_distance_bp'] if nr else None,best_observed_partner_r2=nr['best_r2'] if nr else None,
                      observed_partner_r2_ge_0_8=(nr['best_r2']>=.8 if nr['best_r2'] is not None else None) if nr else None,
                      classification=status,published_variant_identity='NOT_PROVEN_BY_POSITIONAL_MATCH',analysis_scope=SCOPE)
                    output.append(row);counts[(resource,group,policy,status)]+=1
    ctx.check('LD peak denominator',len(peaks),1415);ctx.check('LD peak classification rows three panels',len(output),1415*9)
    ctx.write_csv('ld_peak_coordinate_anchor_classification.csv',output)
    summary=[dict(resource=k[0],sample_group=k[1],pair_policy=k[2],classification=k[3],peak_count=n,peak_denominator=len(peaks),analysis_scope=SCOPE) for k,n in sorted(counts.items())]
    ctx.write_csv('ld_peak_coordinate_anchor_summary.csv',summary)
    return summary

def run(ctx):
    scalar_self_checks(ctx);cfg=ctx.baseline_config['dosage_ld'];lc=cfg['ld']
    bins=ctx.read('dosage_ld/02_ld_distance_profiles.xlsx','distance_bins')
    sums=ctx.read('dosage_ld/05_observed_snp_neighbors.xlsx','summary')
    neigh=ctx.read('dosage_ld/05_observed_snp_neighbors.xlsx','marker_neighbors')
    inv=ctx.read('dosage_ld/02_ld_distance_profiles.xlsx','pair_inventories')
    fb={(r['resource'],r['sample_group'],r['lower_bp_exclusive'],r['upper_bp_inclusive']):r for r in bins if r['coordinate_basis']=='corrected_snp' and r['chromosome']=='ALL'}
    fs={(r['resource'],r['sample_group']):r for r in sums if r['window_bp']==50000}
    fn={(r['resource'],r['sample_group'],r['marker_key']):r for r in neigh if r['window_bp']==50000};fi={r['resource']:r for r in inv}
    decisions=ctx.read('alignment_coverage/01_realignment_evidence.xlsx','marker_decisions');span_index=defaultdict(set)
    for r in decisions:
        if r['accepted_for_mapping']:span_index[str(r['reference_key'])].add((str(r['mapped_chromosome']),int(r['tag_start_bp']),int(r['tag_end_bp']),str(r['mapped_strand'])))
    profiles=[];neighbors_out=[];summaries=[];details=[];tag_audit=[];differences=[];lookup={}
    for resource,panel in ctx.genotype_panels.items():
        print('[LD]',resource,'fresh local pairs and dosage correlations',flush=True)
        pairs=local_pairs(panel.markers,int(lc['max_distance_bp']),int(lc['max_union_pairs']))
        pos=np.array([int(m['position_bp']) for m in panel.markers]);dist=np.abs(pos[pairs[:,1]]-pos[pairs[:,0]])
        fam=[family(m['source_tag_key']) for m in panel.markers]
        spans=[next(iter(span_index[str(m['source_tag_key'])])) if len(span_index[str(m['source_tag_key'])])==1 else None for m in panel.markers]
        same=np.array([fam[i] is not None and fam[i]==fam[j] for i,j in pairs],dtype=bool)
        same_span=np.array([spans[i] is not None and spans[i]==spans[j] for i,j in pairs],dtype=bool)
        ctx.check('LD corrected pair census '+resource,len(pairs),int(fi[resource]['corrected_local_pairs']))
        ctx.check('LD tag family parse '+resource,sum(x is None for x in fam),0)
        ctx.check('LD unique aligned tag span '+resource,sum(x is None for x in spans),0)
        tag_audit.append(dict(resource=resource,qc_markers=len(panel.markers),local_pairs=len(pairs),pairs_same_source_tag_family=int(same.sum()),
          pairs_same_aligned_span=int(same_span.sum()),pairs_same_family_different_span=int(np.sum(same&~same_span)),pairs_same_span_different_family=int(np.sum(same_span&~same)),
          first_bin_pairs=int(np.sum((dist>0)&(dist<=1000))),first_bin_same_family_pairs=int(np.sum(same&(dist>0)&(dist<=1000))),
          family_definition='SOURCE_KEY_PREFIX_BEFORE_TERMINAL_DOUBLE_DASH_SNP_INDEX',corroboration='IDENTICAL_MAPPED_REFERENCE_TAG_START_END_AND_STRAND',analysis_scope=SCOPE))
        values={g:correlation(panel.dosage[:,idx],pairs,cfg['qc']['min_pair_samples'],lc['batch_pairs']) for g,idx in panel.groups.items()}
        edges=[(None,0)]+list(zip(lc['distance_edges_bp'][:-1],lc['distance_edges_bp'][1:]))
        for group,v in values.items():
            for lo,hi in edges:
                base=dist==0 if lo is None else (dist>lo)&(dist<=hi)
                for policy,include in [('ALL_LOCAL_PAIRS',np.ones(len(pairs),bool)),('SAME_SOURCE_TAG_FAMILY_ONLY',same),('EXCLUDE_SAME_SOURCE_TAG_FAMILY',~same),('EXCLUDE_IDENTICAL_ALIGNED_TAG_SPAN',~same_span)]:
                    stats=statistics(v,base&include);profiles.append(dict(resource=resource,sample_group=group,coordinate_basis='corrected_snp',lower_bp_exclusive=lo,upper_bp_inclusive=hi,pair_policy=policy,**stats,analysis_scope=SCOPE))
                    if policy=='ALL_LOCAL_PAIRS':differences.extend(compare(ctx,'LD frozen bin '+resource+':'+group+':'+str(hi),stats,fb[(resource,group,lo,hi)],stats))
            for policy,include in [('ALL_LOCAL_PAIRS',np.ones(len(pairs),bool)),('EXCLUDE_SAME_SOURCE_TAG_FAMILY',~same),('EXCLUDE_IDENTICAL_ALIGNED_TAG_SPAN',~same_span)]:
                rows,summary=neighbors(panel,pairs,dist,v,include)
                fields=list(summary);summary.update(resource=resource,sample_group=group,pair_policy=policy,analysis_scope=SCOPE);summaries.append(summary)
                for i,r in enumerate(rows):lookup[(resource,group,policy,i)]=r;neighbors_out.append(dict(resource=resource,sample_group=group,pair_policy=policy,**r,analysis_scope=SCOPE))
                if policy=='ALL_LOCAL_PAIRS':
                    differences.extend(compare(ctx,'LD frozen50kb summary '+resource+':'+group,summary,fs[(resource,group)],fields));bad=[]
                    for r in rows:
                        ref=fn[(resource,group,r['marker_key'])]
                        for key in ('eligible_other_markers','valid_other_markers','best_partner_key','best_partner_distance_bp','best_r2','status'):
                            if not equal(r[key],ref[key]):bad.append(dict(comparison='50kb marker neighbors '+resource+':'+group,marker_key=r['marker_key'],field=key,observed=r[key],frozen=ref[key]))
                    ctx.check('LD frozen50kb marker rows '+resource+':'+group,len(bad),0);differences.extend(bad)
        for z in np.flatnonzero(((dist>0)&(dist<=1000))|same):
            i,j=map(int,pairs[z]);a,b=panel.markers[i],panel.markers[j]
            row=dict(resource=resource,marker_a=a['marker_key'],marker_b=b['marker_key'],chromosome=str(a['chromosome']),position_a=int(a['position_bp']),position_b=int(b['position_bp']),distance_bp=int(dist[z]),
              source_tag_family_a=fam[i],source_tag_family_b=fam[j],same_source_tag_family=bool(same[z]),same_aligned_tag_span=bool(same_span[z]),first_0_to_1kb_bin=bool(0<dist[z]<=1000),analysis_scope=SCOPE)
            for group,v in values.items():row[group+'_n']=int(v['n'][z]);row[group+'_r2']=scalar(v['r2'][z]);row[group+'_status_code']=int(v['status'][z])
            details.append(row)
    for name,rows in [('ld_same_tag_distance_profiles.csv',profiles),('ld_same_tag_pair_details.csv',details),('ld_tag_family_audit.csv',tag_audit),('ld_observed_neighbors_50kb.csv',neighbors_out),('ld_observed_neighbors_50kb_summary.csv',summaries)]:ctx.write_csv(name,rows)
    ctx.write_csv('ld_frozen_comparison_differences.csv',differences,['comparison','marker_key','field','observed','frozen'])
    peak_tables(ctx,ctx.genotype_panels,lookup)
    ctx.summary['ld_audit']=dict(analysis_scope=SCOPE,same_tag_audit=tag_audit,validation_difference_cells=len(differences),peak_classification_rows=1415*9,
      adjusted_ld_status='DEFERRED_UNTIL_CODING_EVIDENCE_REVIEW',first_bin_profiles=[r for r in profiles if r['upper_bp_inclusive']==1000],neighbor_50kb_summaries=summaries)
    ctx.note('LD is conditional on A/M/C = 0/1/2. Empirical QC does not prove symbol semantics. This module reports unadjusted correlations.')
    ctx.note('A peak coordinate anchor is an observed genotype marker at a published position; variant identity/allele equivalence is unproven. Absent focal genotypes are NOT_ASSESSABLE, never low LD.')
    ctx.note('Focal genotype availability is assessed within each target panel. This module evaluates within-panel focal/partner combinations.')
    ctx.note('Same-source-tag sensitivity removes pairs sharing the identifier prefix before terminal --SNPindex. Identical aligned tag spans provide corroboration and an additional, broader exclusion sensitivity (different source tag identifiers can align to the same span).')
