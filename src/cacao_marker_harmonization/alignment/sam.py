"""Conservative score-separated, biallelic tag placements from Bowtie 2 SAM.

Both allele tags must support the same oriented SNP location. Competition is
assessed BEFORE edit/gap filters, so an inconvenient alternative is not discarded
to manufacture uniqueness. -k is a bounded heuristic search, never an exhaustive
proof of genomic uniqueness. All accepted offsets are derived from CIGAR and
checked against the reference, including reverse-strand cases.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import gzip
import re
from .common import DataError
from .tags import Tag, revcomp

_CIGAR=re.compile(r'(\d+)([MIDNSHP=X])')
@dataclass
class Placement:
    read_name: str
    seqid: str
    start: int
    end: int
    strand: str
    cigar: str
    score: int
    nm: int
    mapq: int
    snp: int|None
    base: str
    other_base: str
    nonfocal_edits: int
    gap_bases: int
    focal_gap_distance: int|None
    reason: str
    reported_xs: int|None=None

    def locus(self):
        return (self.seqid,self.strand,self.snp) if self.snp is not None else (self.seqid,self.strand,self.start,self.end,self.cigar)


def cigar_ops(cigar):
    ops=[(int(n),op) for n,op in _CIGAR.findall(cigar)]
    if not ops or ''.join(str(n)+op for n,op in ops)!=cigar or any(n<=0 for n,_ in ops): raise DataError('Invalid CIGAR: '+cigar)
    return ops


def parse_record(line: str, tag: Tag, allele: str, ref, cfg:dict) -> Placement|None:
    f=line.rstrip('\r\n').split('\t')
    if len(f)<11: raise DataError('Truncated SAM line.')
    try: flag,pos,mapq=int(f[1]),int(f[3]),int(f[4])
    except ValueError as e: raise DataError('Invalid numeric SAM fields.') from e
    if flag & 4:return None
    if flag & (1|2048):raise DataError('Paired or supplementary alignment in an unpaired end-to-end tag run.')
    qseq=tag.seq(allele);strand='-' if flag & 16 else '+'
    if strand=='-':qseq=revcomp(qseq)
    if f[9] not in ('*',qseq):raise DataError('SAM sequence differs from the source tag/orientation.')
    qf=tag.snp0 if strand=='+' else len(qseq)-1-tag.snp0
    opts={}
    for s in f[11:]:
        parts=s.split(':',2)
        if len(parts)==3:
            if parts[0] in opts:raise DataError('Duplicate SAM optional tag.')
            opts[parts[0]]=(parts[1],parts[2])
    for x in ('AS','NM'):
        if x not in opts or opts[x][0]!='i': raise DataError('Missing integer '+x+' in mapped SAM.')
    score,nm=int(opts['AS'][1]),int(opts['NM'][1]);xs=int(opts['XS'][1]) if opts.get('XS',('',))[0]=='i' else None
    ops=cigar_ops(f[5]);qlen=sum(n for n,o in ops if o in 'MIS=X');rlen=sum(n for n,o in ops if o in 'MDN=X')
    if qlen!=len(qseq):raise DataError('CIGAR query length differs from source.')
    rseq=ref.slice(f[2],pos-1,pos-1+rlen).decode('ascii')
    if len(rseq)!=rlen:raise DataError('Alignment falls outside reference sequence.')
    qi=ri=0;mm=0;gaps=0;gap_penalty=0;snp=None;base='';near=None;focal_mm=0;unsupported=False;ambig=False
    for n,op in ops:
        if op in 'M=X':
            a=qseq[qi:qi+n];b=rseq[ri:ri+n]
            mm+=sum(x!=y for x,y in zip(a,b));ambig |= any(x not in 'ACGT' for x in b)
            if qi<=qf<qi+n:
                offset=qf-qi;snp=pos+ri+offset;base=b[offset];focal_mm=int(a[offset]!=base)
            qi+=n;ri+=n
        elif op=='I':
            d=0 if qi<=qf<qi+n else min(abs(qf-qi),abs(qf-(qi+n-1)))
            near=d if near is None else min(near,d);gaps+=n;gap_penalty+=cfg['gap_open']+n*cfg['gap_extend'];qi+=n
        elif op=='D':
            d=min(abs(qf-qi),abs(qf-(qi-1)));near=d if near is None else min(near,d);gaps+=n;gap_penalty+=cfg['gap_open']+n*cfg['gap_extend'];ri+=n
        elif op=='S':qi+=n;unsupported=True
        elif op=='N':ri+=n;unsupported=True
        else:unsupported=True
    if not unsupported and not ambig:
        if nm!=mm+gaps:raise DataError('SAM NM disagrees with reference/CIGAR edit count.')
        if score!=-(mm*cfg['mismatch_penalty']+gap_penalty):raise DataError('SAM AS disagrees with configured end-to-end penalties/reference.')
    alleles=(tag.ref[tag.snp0],tag.alt[tag.snp0]);alleles=tuple(revcomp(x) for x in alleles) if strand=='-' else alleles
    nonfocal=mm+gaps-focal_mm
    reason='ELIGIBLE_PLACEMENT'
    if unsupported:reason='NOT_FULL_END_TO_END'
    elif snp is None:reason='FOCAL_SNP_NOT_ON_REFERENCE_BASE'
    elif ambig or base not in alleles:reason='REFERENCE_ALLELE_UNSUPPORTED'
    elif near is not None and near<=cfg['focal_gap_guard_bases']:reason='GAP_NEAR_FOCAL_SNP'
    elif gaps>cfg['max_gap_bases']:reason='TOO_MANY_GAP_BASES'
    elif nonfocal>cfg['max_nonfocal_edit_bases']:reason='TOO_MANY_NONFOCAL_EDITS'
    other=next((x for x in alleles if x!=base),'')
    return Placement(f[0],f[2],pos,pos+rlen-1,strand,f[5],score,nm,mapq,snp,base,other,nonfocal,gaps,near,reason,xs)


def allele_decision(hits:list[Placement],cfg:dict):
    if not hits:return dict(status='UNMAPPED',reported_hits=0),None
    loci={}
    for h in hits:
        k=h.locus();loci[k]=max(h.score,loci.get(k,-10**9))
    order=sorted(loci.items(),key=lambda x:(-x[1],repr(x[0])))
    best_key,best_score=order[0];next_score=order[1][1] if len(order)>1 else None
    chosen=sorted([h for h in hits if h.locus()==best_key and h.score==best_score],key=lambda h:(h.start,h.end,h.cigar))[0]
    # Also honor XS on a best alignment when the reported list omits that competitor.
    xs_values=[h.reported_xs for h in hits if h.score==best_score and h.reported_xs is not None]
    if xs_values: next_score=max(([next_score] if next_score is not None else [])+xs_values)
    gap=best_score-next_score if next_score is not None else None
    status='RETAINED_BOUNDED_SEARCH'
    if len(hits)>=cfg['max_reported_hits']:status='REPORTING_CAP_REACHED'
    elif gap is not None and gap<cfg['minimum_score_gap']:status='INSUFFICIENT_SCORE_SEPARATION'
    elif best_score-float(cfg['score_min'].split(',')[1])<cfg['minimum_score_gap']:status='SEARCH_THRESHOLD_MARGIN_INSUFFICIENT'
    elif chosen.reason!='ELIGIBLE_PLACEMENT':status=chosen.reason
    elif len({(h.base,h.other_base,h.snp) for h in hits if h.locus()==best_key and h.score==best_score})!=1:status='TIED_CIGAR_FOCAL_CONFLICT'
    row=dict(status=status,reported_hits=len(hits),reported_loci=len(loci),best_score=best_score,next_score=next_score,
             score_gap=gap,best_mapq=chosen.mapq,best_sequence_id=chosen.seqid,best_snp_bp=chosen.snp,
             best_strand=chosen.strand,best_cigar=chosen.cigar,best_NM=chosen.nm,nonfocal_edits=chosen.nonfocal_edits,
             gap_bases=chosen.gap_bases,focal_gap_distance=chosen.focal_gap_distance,
             alternative_scope='Reported end-to-end alignments plus XS; -k/seed search is heuristic, not exhaustive uniqueness.')
    return row,chosen if status=='RETAINED_BOUNDED_SEARCH' else None


def combine_alleles(tag:Tag,dr,hr,da,ha,chromosomes):
    row=dict(query_id=tag.query_id,source_row=tag.source_row,reference_key=tag.reference_key,source_marker_id=tag.marker_id,
             ref_allele_status=dr['status'],alt_allele_status=da['status'],accepted_for_mapping=False,
             placement_status='BOTH_ALLELES_REQUIRED')
    if hr is None or ha is None:return row
    if hr.locus()!=ha.locus():row['placement_status']='ALLELE_PLACEMENTS_DISAGREE';return row
    if hr.base!=ha.base or hr.other_base!=ha.other_base:raise DataError('Paired allele reference bases disagree.')
    chrom=chromosomes.get(hr.seqid)
    row.update(placement_status='RETAINED_PRIMARY_CHROMOSOME' if chrom else 'RETAINED_UNPLACED_CONTIG',
        accepted_for_mapping=bool(chrom),mapped_sequence_id=hr.seqid,mapped_chromosome=chrom,mapped_snp_bp=hr.snp,
        mapped_strand=hr.strand,tag_start_bp=hr.start,tag_end_bp=hr.end,alt_tag_start_bp=ha.start,alt_tag_end_bp=ha.end,
        ref_cigar=hr.cigar,alt_cigar=ha.cigar,ref_score=hr.score,alt_score=ha.score,
        ref_nonfocal_edits=hr.nonfocal_edits,alt_nonfocal_edits=ha.nonfocal_edits,
        assembly_base=hr.base,other_source_base=hr.other_base)
    return row


def parse_sam(path,tags,ref,cfg,chromosomes,log=None):
    results={};alternates=[];seen=set();sq={};name=None;lines=[];records=0
    def finish(n,ls):
        if n is None:return
        if n in seen:raise DataError('Noncontiguous/repeated QNAME in SAM; use --reorder.')
        seen.add(n)
        if '__' not in n:raise DataError('Unknown SAM QNAME.')
        qid,allele=n.rsplit('__',1)
        if qid not in tags or allele not in ('R','A'):raise DataError('Unknown SAM tag ID.')
        hits=[h for s in ls if (h:=parse_record(s,tags[qid],allele,ref,cfg)) is not None]
        if len(ls)>cfg['max_reported_hits']:raise DataError('SAM exceeds configured reporting cap.')
        if any(int(s.split('\t')[1])&4 for s in ls) and len(ls)!=1:raise DataError('Mixed mapped/unmapped SAM group.')
        dec,best=allele_decision(hits,cfg);results[n]=(dict(query_id=qid,allele=allele,**dec),best)
        # Every alignment remains in the cached compressed SAM. The workbook shows
        # the two highest-score records per allele, not a falsely complete hit list.
        for rank,h in enumerate(sorted(hits,key=lambda h:(-h.score,h.seqid,h.start,h.cigar))[:2],1):
            alternates.append(dict(query_id=qid,allele=allele,reported_rank=rank,**asdict(h)))
    op=gzip.open if path.suffix=='.gz' else open
    with op(path,'rt',encoding='ascii') as f:
        for line in f:
            if line.startswith('@'):
                if line.startswith('@SQ\t'):
                    kv=dict(x.split(':',1) for x in line.strip().split('\t')[1:])
                    if kv['SN'] in sq:raise DataError('Duplicate SAM reference sequence ID.')
                    sq[kv['SN']]=int(kv['LN'])
                continue
            n=line.split('\t',1)[0]
            if name!=n:finish(name,lines);name=n;lines=[]
            lines.append(line);records+=1
            if len(lines)>cfg['max_reported_hits']:raise DataError('Excessive alignments in one SAM group.')
    finish(name,lines)
    if sq!={k:v['length'] for k,v in ref.entries.items()}:raise DataError('SAM reference sequence dictionary differs from validated genome.')
    expected={q+'__'+a for q in tags for a in ('R','A')}
    if seen!=expected:raise DataError('SAM missing or unexpected allele records.')
    combined=[]
    for qid,tag in tags.items():
        dr,hr=results[qid+'__R'];da,ha=results[qid+'__A'];combined.append(combine_alleles(tag,dr,hr,da,ha,chromosomes))
    return combined,[v[0] for v in results.values()],alternates,records
