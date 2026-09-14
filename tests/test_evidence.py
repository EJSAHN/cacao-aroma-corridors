import unittest
import tempfile
from pathlib import Path
from collections import Counter
from cacao_marker_harmonization.evidence.common import *
from cacao_marker_harmonization.evidence.sequences import *
from cacao_marker_harmonization.evidence.events import *
from cacao_marker_harmonization.evidence.genotypes import *
from cacao_marker_harmonization.evidence.readers import Matrix


def base_event(**kw):
    r=dict(source_file='source.xlsx',source_sheet='master',source_row=3,chromosome='1',peak_bp=100,trait='trait_R',
      method_reported='GLM',marker_filter_reported='G7',source_start_bp=90,source_end_bp=110,source_block_number=2,
      p_value_reported=.0001,explained_variance_reported=.15,event_key='key',in_master=True,raw_record='{}')
    r.update(kw);return r

def event_cfg(**kw):
    r=dict(chromosomes=['1'],association_master_sheet='master',trait_alias_candidates=[],reported_variance_scale_candidates=[1,100,.01])
    r.update(kw);return r

def seqrow(**kw):
    r=dict(source_row=2,marker_id='1|F|0--2',marker_key='1|F|0--2',reference_key='1|F|0--2',reference_sequence_name='1|F|0--2',
        chromosome='1',position_bp=100,five_flank='AC',three_flank='GT',variation='G>T',assembly_reported='V2',remark='',strand_reported='')
    r.update(kw);return r

def reference(s='ACGGT',row=2):return dict(sequence=s,source_row=row)

class ValueTests(unittest.TestCase):
    def test_bare_accession(self):self.assertEqual(chromosome('NC_030859.1'),'NC_030859.1')
    def test_explicit_chromosome(self):self.assertEqual(chromosome('NC_030859.1_chromosome_10'),'10')
    def test_dart_spelling(self):self.assertEqual(canonical_marker('123.F.0.5'),'123|F|0--5')
    def test_no_generic_punctuation_removal(self):self.assertNotEqual(canonical_marker('ab.cd'),canonical_marker('ab_cd'))
    def test_integer(self):self.assertEqual(integer('1 234'),1234);self.assertIsNone(integer('1.25'));self.assertIsNone(integer('#VALUE!'))
    def test_zero_p_not_missing(self):self.assertEqual(number('0'),0)
    def test_suffix(self):self.assertEqual(dart_suffix('123.F.0.5'),5);self.assertIsNone(dart_suffix('mk5'))

class SequenceTests(unittest.TestCase):
    def test_reconstruct_reference(self):self.assertEqual(reconstruct(seqrow())['ref_tag'],'ACGGT')
    def test_reconstruct_alternate(self):self.assertEqual(reconstruct(seqrow())['alt_tag'],'ACTGT')
    def test_empty_flank_allowed(self):self.assertEqual(reconstruct(seqrow(five_flank=''))['ref_tag'],'GGT')
    def test_no_sequence_not_accepted(self):self.assertIsNone(reconstruct(seqrow(five_flank='',three_flank='')))
    def test_indel_rejected(self):self.assertIsNone(reconstruct(seqrow(variation='A>AT')))
    def test_ambiguous_flank_rejected(self):self.assertIsNone(reconstruct(seqrow(five_flank='ANC')))
    def test_variant_direction_retained(self):self.assertIsNone(reconstruct(seqrow(variation='A/A')))
    def test_identical_alleles_rejected(self):self.assertIsNone(reconstruct(seqrow(variation='A>A')))
    def test_exact(self):self.assertEqual(compare_tag(reconstruct(seqrow()),[reference()]),'exact_ref_tag')
    def test_alternate_distinct(self):self.assertEqual(compare_tag(reconstruct(seqrow()),[reference('ACTGT')]),'exact_alt_tag')
    def test_reverse(self):self.assertEqual(compare_tag(reconstruct(seqrow()),[reference('ACCGT')]),'reverse_complement_ref_tag')
    def test_mismatch(self):self.assertEqual(compare_tag(reconstruct(seqrow()),[reference('AAAAA')]),'sequence_mismatch')
    def test_missing(self):self.assertEqual(compare_tag(reconstruct(seqrow()),[]),'missing_reference_name')
    def test_duplicate_same_sequence_not_silently_used(self):self.assertEqual(compare_tag(reconstruct(seqrow()),[reference(),reference(row=3)]),'duplicate_reference_name')
    def test_two_keys_compared(self):
        r=seqrow(marker_key='2|F|0--2',marker_id='2|F|0--2')
        d,s,a=sequence_evidence([r],{'1|F|0--2':[reference()],'2|F|0--2':[reference('AAAAA')]})
        self.assertTrue(d[0]['sequence_supported_reference_key']);self.assertEqual(d[0]['marker_lookup_status'],'sequence_mismatch')
    def test_models_are_separate(self):
        d,_,_=sequence_evidence([seqrow()],{'1|F|0--2':[reference()]});m=models(d[0])
        self.assertEqual(m['five_flank_start_1based'],98);self.assertEqual(m['five_flank_start_0based'],97)
        self.assertEqual(m['suffix_forward_0'],98);self.assertEqual(m['suffix_reverse_0'],97)
    def test_unplaced_coordinates(self):self.assertEqual(coordinate_relation(dict(chromosome='ctg',position_bp=3),dict(chromosome='1',position_bp=3),{'1'}),'unavailable')
    def test_different_chrom_not_offset(self):self.assertEqual(coordinate_relation(dict(chromosome='2',position_bp=3),dict(chromosome='1',position_bp=3),{'1','2'}),'different_chromosomes')

class EventTests(unittest.TestCase):
    def test_suffix_whitespace(self):self.assertEqual(lexical_trait('acetate _R'),'acetate_R')
    def test_ur_not_r(self):self.assertNotEqual(lexical_trait('acetate_UR'),lexical_trait('acetate_R'))
    def test_internal_identity_preserved(self):self.assertNotEqual(lexical_trait('a  b_R'),lexical_trait('a b_R'))
    def test_different_model_preserved(self):self.assertNotEqual(event_tuple(base_event()),event_tuple(base_event(method_reported='MLM')))
    def test_extra_alias_only_candidate(self):
        rows=[base_event(),base_event(source_sheet='extra',in_master=False,trait='trait _R',source_row=4)]
        r=reconcile(rows,event_cfg());self.assertEqual(r['crosswalk'][1]['matching_basis'],'terminal_suffix_whitespace_candidate')
        self.assertEqual(rows[1]['trait'],'trait _R');self.assertEqual(len(r['crosswalk']),2)
    def test_p_difference_exposed(self):
        r=reconcile([base_event(),base_event(source_sheet='extra',in_master=False,p_value_reported=.001)],event_cfg())
        self.assertEqual(r['crosswalk'][1]['p_comparison'],'different');self.assertEqual(r['crosswalk'][1]['action'],'NO_PRIMARY_CHANGE')
    def test_percent_difference_not_applied(self):
        rr=base_event(source_sheet='extra',in_master=False,explained_variance_reported=15)
        r=reconcile([base_event(),rr],event_cfg());self.assertEqual(r['crosswalk'][1]['source_to_master_variance_factor'],100);self.assertEqual(rr['explained_variance_reported'],15)
    def test_rotated_interval_reported(self):
        r=reconcile([base_event(),base_event(source_sheet='extra',in_master=False,source_start_bp=2,source_end_bp=90,source_block_number=110)],event_cfg())
        self.assertEqual(r['crosswalk'][1]['interval_relation'],'field_permutation_matches_master')
    def test_missing_boundary_not_imputed(self):
        rr=base_event(source_start_bp=None,source_end_bp=None)
        r=reconcile([rr],event_cfg());self.assertEqual(r['interval_review'][0]['state'],'no_boundaries');self.assertIsNone(rr['source_start_bp'])
    def test_ambiguous_master_not_selected(self):
        r=reconcile([base_event(),base_event(source_row=4),base_event(source_sheet='extra',in_master=False)],event_cfg())
        self.assertEqual(r['crosswalk'][2]['match_status'],'ambiguous_master_event');self.assertIsNone(r['crosswalk'][2]['master_row'])
    def test_alias_cannot_apply(self):
        with self.assertRaises(DataError):validated_aliases({'trait_alias_candidates':[dict(source='a',target='b',reason='test',apply_to_primary=True)]})
    def test_alias_needs_reason(self):
        with self.assertRaises(DataError):validated_aliases({'trait_alias_candidates':[dict(source='a',target='b')]})
    def test_alias_chain_rejected(self):
        with self.assertRaises(DataError):validated_aliases({'trait_alias_candidates':[dict(source='a',target='b',reason='x'),dict(source='b',target='c',reason='y')]})
    def test_config_alias_p_payload_separate(self):
        cfg=event_cfg(trait_alias_candidates=[dict(source='other',target='trait_R',reason='review',apply_to_primary=False)])
        r=reconcile([base_event(),base_event(source_sheet='extra',in_master=False,trait='other',p_value_reported=.002)],cfg)
        self.assertEqual(r['crosswalk'][1]['matching_basis'],'configured_trait_alias_candidate');self.assertEqual(r['crosswalk'][1]['p_comparison'],'different')
    def test_tiny_p_not_rounded_to_equal(self):self.assertEqual(scalar_comparison(1e-9,1e-8,[1],1e-9)[0],'different')

class GenotypeTests(unittest.TestCase):
    def test_GTM_conflict(self):
        a,b,c=interpret_symbols(Counter('GTMM'),{'N',''})
        self.assertEqual(a,'more_than_two_literal_alleles');self.assertEqual(b,'two_base_labels_plus_M')
    def test_ACM_literal(self):self.assertEqual(interpret_symbols(Counter('ACMM'),{'N',''})[0],'compatible_with_literal_iupac')
    def test_one_base_M_is_unresolved(self):self.assertEqual(interpret_symbols(Counter('AMM'),{'N'})[1],'M_meaning_unidentifiable_from_observed_labels')
    def test_unknown_token(self):self.assertEqual(interpret_symbols(Counter({'X':3}),{'N'})[0],'unsupported_symbols')
    def test_missing_not_allele(self):self.assertEqual(interpret_symbols(Counter({'N':3,'':2}),{'N',''})[0],'no_calls')
    def test_more_than_two_bases(self):self.assertEqual(interpret_symbols(Counter('ACGM'),{'N'})[1],'more_than_two_observed_base_labels')
    def setup_pair(self,calls_a=('A','M','N'),calls_b=('A','C','N'),chrom_b='2'):
        ma=dict(source_row=2,canonical_id='x',marker_id='x',chromosome='1',position_bp=9,record_index=0,alleles_reported='A>C')
        mb={**ma,'chromosome':chrom_b}
        a=Matrix('a',['s1','s2','s3'],[ma],[calls_a],'s','reader',[])
        b=Matrix('b',['s1','s2','s3'],[mb],[calls_b],'s','reader',[])
        cfg=dict(chromosomes=['1','2'],missing_symbols=['N',''],symbol_polarity_swap={'A':'C','C':'A','M':'M'})
        return a,b,cfg
    def test_concordance_not_location_validation(self):
        a,b,cfg=self.setup_pair();r=compare_matrices(a,b,cfg)
        d=r['matched_marker_calls'][0];self.assertEqual(d['raw_symbol_concordance'],.5);self.assertEqual(len(r['chromosome_conflicts']),1)
    def test_swap_separate(self):
        a,b,cfg=self.setup_pair(('A','M','C'),('C','M','A'));r=compare_matrices(a,b,cfg)['matched_marker_calls'][0]
        self.assertAlmostEqual(r['raw_symbol_concordance'],1/3);self.assertEqual(r['swapped_concordance'],1)
        self.assertEqual(a.calls[0],('A','M','C'))
    def test_no_pairs_fraction_missing(self):
        a,b,cfg=self.setup_pair(('N','N','N'),('N','N','N'));self.assertIsNone(compare_matrices(a,b,cfg)['matched_marker_calls'][0]['raw_symbol_concordance'])
    def test_duplicate_marker_not_selected(self):
        a,b,cfg=self.setup_pair();a.markers.append(dict(a.markers[0]));r=compare_matrices(a,b,cfg)
        self.assertEqual(len(r['matched_marker_calls']),0);self.assertEqual(len(r['ambiguous_marker_ids']),1)

if __name__=='__main__':unittest.main()
