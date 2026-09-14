"""Read the immutable, sequence-backed source tags; do not reinterpret genotypes."""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import re
from .common import DataError, integer, text

_COMPLEMENT = str.maketrans('ACGTN','TGCAN')
def revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]

@dataclass(frozen=True)
class Tag:
    query_id: str
    source_row: int
    reference_key: str
    marker_id: str
    ref: str
    alt: str
    snp0: int

    def seq(self, allele: str):
        if allele not in ('R','A'): raise DataError('Invalid source allele selector.')
        return self.ref if allele == 'R' else self.alt


def source_tags(evidence: list[dict], prior: list[dict], min_length=25):
    old={r['query_id']:r for r in prior}
    if len(old)!=len(prior): raise DataError('Duplicate prior tag query IDs.')
    tags={};audit=[];seen=set()
    for r in evidence:
        rn=integer(r.get('source_row'))
        if rn is None or rn<1 or rn in seen: raise DataError('Invalid/duplicate tag source row.')
        seen.add(rn); qid='diversity_row_'+str(rn); p=old.get(qid)
        if p is None or text(p['source_marker_id'])!=text(r['marker_id']) or text(p['reference_key'])!=text(r['reference_key']):
            raise DataError('Source tag inventory differs from exact-mapping inventory.')
        state=text(p['eligibility'])
        out={'query_id':qid,'source_row':rn,'marker_id':text(r['marker_id']),'reference_key':text(r['reference_key']),
             'prior_mapping_status':p['mapping_status'],'input_status':state}
        if state=='ELIGIBLE':
            left=re.sub(r'\s+','',text(r['five_flank'])).upper();right=re.sub(r'\s+','',text(r['three_flank'])).upper()
            var=re.fullmatch(r'([ACGT])\s*>\s*([ACGT])',text(r['variation']).upper())
            if not var or var[1]==var[2] or not re.fullmatch('[ACGT]*',left+right): raise DataError('Prior-eligible tag is invalid.')
            ref=left+var[1]+right;alt=left+var[2]+right
            if ref!=re.sub(r'\s+','',text(r['reference_sequence'])).upper() or r['ref_lookup_status']!='exact_ref_tag':
                raise DataError('Prior-eligible source tag differs from its reference-sequence evidence.')
            if hashlib.sha256(ref.encode()).hexdigest()!=p['sequence_sha256'] or len(left)!=integer(p['source_snp_index_0based']) or len(ref)!=integer(p['tag_length']):
                raise DataError('Tag hash/length/SNP offset differs from the checked exact-tag inventory.')
            if len(ref)<min_length: out['input_status']='BELOW_CONFIGURED_MINIMUM_LENGTH'
            else: tags[qid]=Tag(qid,rn,text(r['reference_key']),text(r['marker_id']),ref,alt,len(left))
        audit.append(out)
    if set(old)!=set(r['query_id'] for r in audit): raise DataError('Incomplete sequence evidence.')
    if not tags: raise DataError('No eligible sequence tags.')
    return tags,audit


def write_fasta(path,tags):
    with path.open('x',encoding='ascii',newline='\n') as f:
        for q in tags.values():
            for allele in ('R','A'): f.write(f'>{q.query_id}__{allele}\n{q.seq(allele)}\n')
