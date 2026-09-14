"""Symbol diagnostics and matched-sample concordance without assigning causal alleles."""
from __future__ import annotations
from collections import Counter
from .common import fraction
from .sequences import make_index,coordinate_relation

IUPAC={'A':{'A'},'C':{'C'},'G':{'G'},'T':{'T'},'M':{'A','C'},'R':{'A','G'},
       'W':{'A','T'},'S':{'C','G'},'Y':{'C','T'},'K':{'G','T'}}


def interpret_symbols(counts,missing):
    observed=set(counts)-set(missing)
    unsupported=observed-set(IUPAC)
    alleles=set().union(*(IUPAC.get(x,set()) for x in observed))
    literal='unsupported_symbols' if unsupported else 'no_calls' if not observed else 'more_than_two_literal_alleles' if len(alleles)>2 else 'compatible_with_literal_iupac'
    base=observed&set('ACGT');other=observed-(set('ACGT')|{'M'})
    generic='unsupported_symbols' if other else 'more_than_two_observed_base_labels' if len(base)>2 else 'no_calls' if not observed else 'two_base_labels_plus_M' if len(base)==2 and 'M' in observed else 'M_meaning_unidentifiable_from_observed_labels' if 'M' in observed else 'base_labels_only'
    return literal,generic,''.join(sorted(alleles))


def symbol_diagnostics(matrices,cfg):
    marker=[];resource=[];sample=[];symbols=[]
    missing=set(cfg['missing_symbols'])
    for name,matrix in matrices.items():
        totals=Counter();patterns=Counter();samplecounts=[Counter() for _ in matrix.samples]
        for record,calls in zip(matrix.markers,matrix.calls):
            co=Counter(calls);totals.update(co)
            for cc,t in zip(samplecounts,calls):cc[t]+=1
            literal,generic,alleles=interpret_symbols(co,missing);patterns[(literal,generic)]+=1
            marker.append(dict(resource=name,source_row=record['source_row'],marker_id=record['marker_id'],
                observed_symbols=';'.join(sorted(co)),observed_counts='; '.join(f'{k or "<blank>"}:{n}' for k,n in sorted(co.items())),
                missing_calls=sum(n for k,n in co.items() if k in missing),M_calls=co['M'],total_calls=len(calls),
                literal_iupac_state=literal,literal_alleles=alleles,generic_M_hypothesis=generic,
                biological_encoding_confirmed=False))
        for (literal,generic),n in sorted(patterns.items()):
            resource.append(dict(resource=name,literal_iupac_state=literal,generic_M_hypothesis=generic,markers=n,
                evidence_status='SYMBOL_COUNTS_ONLY_NOT_A_CODEBOOK'))
        symbols.extend(dict(resource=name,symbol=k or '<blank>',calls=n,is_missing=k in missing) for k,n in sorted(totals.items()))
        for sid,co in zip(matrix.samples,samplecounts):
            n=sum(co.values());nm=sum(v for k,v in co.items() if k in missing)
            sample.append(dict(resource=name,sample_id=sid,marker_rows=n,missing_calls=nm,non_missing_call_fraction=fraction(n-nm,n),
                M_calls=co['M'],M_fraction_of_nonmissing=fraction(co['M'],n-nm)))
    return dict(marker_symbols=marker,resource_patterns=resource,sample_symbols=sample,symbol_counts=symbols)


def compare_matrices(a,b,cfg):
    ai=make_index(a.markers,'canonical_id');bi=make_index(b.markers,'canonical_id')
    a_sample={s:i for i,s in enumerate(a.samples)};b_sample={s:i for i,s in enumerate(b.samples)}
    common=sorted(set(a.samples)&set(b.samples));pairs=[(a_sample[s],b_sample[s]) for s in common]
    missing=set(cfg['missing_symbols']);swap=cfg['symbol_polarity_swap'];chromosomes=set(cfg['chromosomes'])
    summaries=Counter();detail=[];discord=[];ambiguities=[]
    for k in sorted(set(ai)&set(bi)):
        if len(ai[k])!=1 or len(bi[k])!=1:
            ambiguities.append(dict(canonical_id=k,rows_a=len(ai[k]),rows_b=len(bi[k]),action='NOT_COMPARED'));continue
        ra=ai[k][0];rb=bi[k][0];ca=a.calls[ra['record_index']];cb=b.calls[rb['record_index']]
        relation=coordinate_relation(ra,rb,chromosomes)
        n=equal=swapequal=swapeligible=bothmiss=onemiss=0
        for sid,(ia,ib) in zip(common,pairs):
            x=ca[ia];y=cb[ib]
            if x in missing or y in missing:
                bothmiss+=int(x in missing and y in missing);onemiss+=int((x in missing)!=(y in missing));continue
            n+=1;equal+=int(x==y)
            if x in swap and y in swap:
                swapeligible+=1;swapequal+=int(x==swap[y])
            if x!=y:
                discord.append(dict(canonical_id=k,sample_id=sid,symbol_a=x,symbol_b=y,coordinate_relation=relation,
                    source_row_a=ra['source_row'],source_row_b=rb['source_row']))
        row=dict(resource_a=a.name,resource_b=b.name,canonical_id=k,source_row_a=ra['source_row'],source_row_b=rb['source_row'],
            chromosome_a=ra['chromosome'],position_a_bp=ra['position_bp'],chromosome_b=rb['chromosome'],position_b_bp=rb['position_bp'],
            coordinate_relation=relation,common_samples=len(common),paired_nonmissing=n,equal_symbols=equal,different_symbols=n-equal,
            raw_symbol_concordance=fraction(equal,n),both_missing=bothmiss,one_missing=onemiss,
            swap_eligible_pairs=swapeligible,equal_after_configured_swap=swapequal,swapped_concordance=fraction(swapequal,swapeligible),
            source_alleles_a=ra['alleles_reported'],source_alleles_b=rb['alleles_reported'])
        detail.append(row)
        summaries[(relation,'markers')]+=1;summaries[(relation,'paired')]+=n;summaries[(relation,'equal')]+=equal
    summary=[dict(resource_a=a.name,resource_b=b.name,coordinate_relation=rel,markers=summaries[(rel,'markers')],
                  paired_nonmissing=summaries[(rel,'paired')],equal_symbols=summaries[(rel,'equal')],
                  raw_symbol_concordance=fraction(summaries[(rel,'equal')],summaries[(rel,'paired')])) for rel in sorted({r['coordinate_relation'] for r in detail})]
    sample_overlap=[dict(sample_id=s,in_a=s in a_sample,in_b=s in b_sample) for s in sorted(set(a.samples)|set(b.samples))]
    return dict(pair_summary=summary,matched_marker_calls=detail,chromosome_conflicts=[r for r in detail if r['coordinate_relation']=='different_chromosomes'],
                discordant_cells=discord,ambiguous_marker_ids=ambiguities,sample_overlap=sample_overlap)
