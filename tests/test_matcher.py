import random
import unittest
from cacao_marker_harmonization.reference.common import DataError
from cacao_marker_harmonization.reference.matcher import (Query, Hit, ExactTagMatcher, reverse_complement,
    encode_seed, select_seed, build_queries, summarize_mappings)


def make_query(ref=b'GACTTACGAGTCCTAGCGATCACTG', ix=12, qid='q'):
    other = next(b for b in b'ACGT' if b != ref[ix])
    return Query(qid, 2, 'key', ref, ref[:ix] + bytes([other]) + ref[ix + 1:], ix)


def brute(queries, seqid, seq):
    out=[]
    for q in queries:
        for strand in ('+', '-'):
            r=q.ref_tag if strand=='+' else reverse_complement(q.ref_tag)
            a=q.alt_tag if strand=='+' else reverse_complement(q.alt_tag)
            ix=q.snp_index if strand=='+' else len(r)-1-q.snp_index
            for start in range(len(seq)-len(r)+1):
                segment=seq[start:start+len(r)]
                if segment==r:
                    tag='source_ref';base=r[ix];other=a[ix]
                elif segment==a:
                    tag='source_alt';base=a[ix];other=r[ix]
                else:continue
                out.append(Hit(q.query_id,seqid,start+1,start+len(r),start+ix+1,strand,tag,chr(base),chr(other)))
    return set(out)


class ExactMatchingTests(unittest.TestCase):
    def scan(self, q, seq, **kwargs):
        return list(ExactTagMatcher([q],seed_k=4,chunk_bases=100,**kwargs).scan('I',seq))
    def test_forward_reference(self):
        q=make_query();h=self.scan(q,b'NN'+q.ref_tag+b'NN')
        self.assertEqual(h,[Hit('q','I',3,27,15,'+','source_ref','C','A')])
    def test_forward_alternate(self):
        q=make_query();h=self.scan(q,q.alt_tag)
        self.assertEqual(len(h),1);self.assertEqual(h[0].matched_source_allele,'source_alt')
        self.assertEqual(h[0].assembly_base,chr(q.alt_tag[q.snp_index]))
    def test_reverse_reference_coordinate(self):
        q=make_query();h=self.scan(q,b'NNNNN'+reverse_complement(q.ref_tag))
        self.assertEqual(len(h),1);self.assertEqual(h[0].strand,'-')
        self.assertEqual(h[0].snp_position_bp,5+len(q.ref_tag)-q.snp_index)
        self.assertEqual(h[0].assembly_base,reverse_complement(q.ref_tag[q.snp_index:q.snp_index+1]).decode())
    def test_reverse_alternate_coordinate(self):
        q=make_query(ix=1);h=self.scan(q,b'N'+reverse_complement(q.alt_tag))
        self.assertEqual(h[0].snp_position_bp,len(q.ref_tag));self.assertEqual(h[0].matched_source_allele,'source_alt')
    def test_index_zero(self):
        q=make_query(ix=0);h=self.scan(q,q.ref_tag)
        self.assertEqual(h[0].snp_position_bp,1)
    def test_last_index_reverse(self):
        q=make_query(ix=24);h=self.scan(q,reverse_complement(q.ref_tag))
        self.assertEqual(h[0].snp_position_bp,1)
    def test_both_alleles_separate_loci(self):
        q=make_query();h=self.scan(q,q.ref_tag+b'NNNN'+q.alt_tag)
        self.assertEqual(len(h),2)
    def test_multimapping_is_not_accepted(self):
        q=make_query();h=self.scan(q,q.ref_tag+b'NNNN'+q.ref_tag)
        aud=[dict(query_id='q',eligibility='ELIGIBLE',source_chromosome='1',source_position_bp=13)]
        r=summarize_mappings(aud,h,{'I':'1'})[0]
        self.assertFalse(r['unique_exact_placement']);self.assertEqual(r['exact_placements'],2)
    def test_unplaced_second_hit_prevents_unique(self):
        q=make_query();m=ExactTagMatcher([q],seed_k=4,chunk_bases=100)
        h=list(m.scan('I',q.ref_tag))+list(m.scan('unplaced',q.ref_tag))
        r=summarize_mappings([dict(query_id='q',eligibility='ELIGIBLE')],h,{'I':'1'})[0]
        self.assertEqual(r['mapping_status'],'MULTIPLE_EXACT_FULL_TAG_PLACEMENTS')
    def test_chunk_boundary_once(self):
        q=make_query();seq=b'N'*92+q.ref_tag+b'N'*80
        self.assertEqual(set(self.scan(q,seq)),brute([q],'I',seq));self.assertEqual(len(self.scan(q,seq)),1)
    def test_seed_at_final_valid_position(self):
        q=make_query(ix=0);seq=b'NCGT'*70+q.ref_tag
        self.assertEqual(set(self.scan(q,seq)),brute([q],'I',seq))
    def test_n_outside_seed_rejects_full_tag(self):
        q=make_query();off,_=select_seed(q.ref_tag,q.snp_index,4)
        ix=next(i for i in range(len(q.ref_tag)) if not off<=i<off+4)
        seq=q.ref_tag[:ix]+b'N'+q.ref_tag[ix+1:]
        self.assertEqual(self.scan(q,seq),[])
    def test_mismatch_not_at_snp_rejected(self):
        q=make_query();seq=b'T'+q.ref_tag[1:]
        self.assertEqual(self.scan(q,seq),[])
    def test_third_allele_rejected(self):
        q=make_query();b=next(x for x in b'ACGT' if x not in (q.ref_tag[q.snp_index],q.alt_tag[q.snp_index]))
        seq=q.ref_tag[:q.snp_index]+bytes([b])+q.ref_tag[q.snp_index+1:]
        self.assertEqual(self.scan(q,seq),[])
    def test_lowercase_reference(self):
        q=make_query();self.assertEqual(len(self.scan(q,q.ref_tag.lower())),1)
    def test_short_reference(self):
        self.assertEqual(self.scan(make_query(),b'AC'),[])
    def test_hit_limit_fails_not_truncates(self):
        q=make_query()
        with self.assertRaises(DataError):self.scan(q,q.ref_tag+b'N'+q.ref_tag,max_hits=1)
    def test_duplicate_query_ids_rejected(self):
        with self.assertRaises(DataError):ExactTagMatcher([make_query(),make_query()],seed_k=4)
    def test_extra_variant_query_rejected(self):
        q=make_query();bad=Query('q',2,'k',q.ref_tag,b'T'+q.alt_tag[1:],q.snp_index)
        with self.assertRaises(DataError):ExactTagMatcher([bad],seed_k=4)
    def test_palindrome_multiple_orientations_not_forced(self):
        q=make_query(ref=b'ACGTACGTACGTACGTACGTACGT',ix=11)
        h=self.scan(q,q.ref_tag)
        self.assertEqual({x.strand for x in h},{'+','-'})
        r=summarize_mappings([dict(query_id='q',eligibility='ELIGIBLE')],h,{'I':'1'})[0]
        self.assertFalse(r['unique_exact_placement'])
    def test_k_bound(self):
        for k in (3,13):
            with self.assertRaises(DataError):ExactTagMatcher([make_query()],seed_k=k)
    def test_seed_excludes_snp(self):
        q=make_query();off,seed=select_seed(q.ref_tag,q.snp_index,4)
        self.assertFalse(off<=q.snp_index<off+4);self.assertEqual(seed,q.ref_tag[off:off+4])
    def test_seed_encoding(self):
        self.assertEqual(encode_seed(b'ACGT'),27)
    def test_bad_seed(self):
        with self.assertRaises(DataError):encode_seed(b'ACNT')
    def test_randomized_against_independent_bruteforce(self):
        rng=random.Random(714)
        for rep in range(40):
            seq=bytearray(rng.choice(b'ACGTN') for _ in range(700))
            qs=[]
            for i in range(8):
                ref=bytes(rng.choice(b'ACGT') for _ in range(rng.randrange(15,41)))
                q=make_query(ref,rng.randrange(len(ref)),f'q{i}')
                qs.append(q)
                ins=q.ref_tag if i%2==0 else q.alt_tag
                if i%3==0:ins=reverse_complement(ins)
                at=20+i*80;seq[at:at+len(ins)]=ins
            seq=bytes(seq)
            got=set(ExactTagMatcher(qs,seed_k=4,chunk_bases=100).scan('x',seq))
            self.assertEqual(got,brute(qs,'x',seq),(rep,len(got)))


class SourceQueryTests(unittest.TestCase):
    def row(self):
        q=make_query()
        return dict(source_row=2,reference_key='a',marker_id='b',chromosome='1',position_bp=3,
                    five_flank=q.ref_tag[:q.snp_index].decode(),three_flank=q.ref_tag[q.snp_index+1:].decode(),
                    variation='C>A',reference_sequence=q.ref_tag.decode(),ref_lookup_status='exact_ref_tag')
    def test_source_query(self):
        qs,au=build_queries([self.row()],25,4);self.assertEqual(len(qs),1);self.assertEqual(au[0]['eligibility'],'ELIGIBLE')
    def test_mismatch_not_assumed_correct(self):
        r=self.row();r['reference_sequence']='T'+r['reference_sequence'][1:]
        q,a=build_queries([r],25,4);self.assertFalse(q);self.assertEqual(a[0]['eligibility'],'FLANK_REFERENCE_TAG_MISMATCH')
    def test_ambiguous_base_not_converted(self):
        r=self.row();r['five_flank']='N'+r['five_flank'][1:]
        q,a=build_queries([r],25,4);self.assertFalse(q)
    def test_prior_status_required(self):
        r=self.row();r['ref_lookup_status']='sequence_mismatch';self.assertFalse(build_queries([r],25,4)[0])
    def test_duplicate_row_rejected(self):
        with self.assertRaises(DataError):build_queries([self.row(),self.row()],25,4)
    def test_minimum_length(self):
        self.assertEqual(build_queries([self.row()],26,4)[1][0]['eligibility'],'TAG_TOO_SHORT')
    def test_unplaced_unique_preserved(self):
        q=make_query();h=list(ExactTagMatcher([q],seed_k=4).scan('unknown_scaffold',q.ref_tag))
        row=dict(query_id='q',eligibility='ELIGIBLE',source_chromosome='1',source_position_bp=1)
        m=summarize_mappings([row],h,{'I':'1'})[0]
        self.assertTrue(m['unique_exact_placement']);self.assertIsNone(m['mapped_chromosome'])
