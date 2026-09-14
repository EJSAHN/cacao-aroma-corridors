import unittest
from cacao_marker_harmonization.reference.common import DataError
from cacao_marker_harmonization.reference.harmonize import project_markers,candidate_checks,peak_links

class HarmonizationTests(unittest.TestCase):
    def mapping(self):
        return dict(source_row=2,reference_key='k1',source_marker_id='wrong_name',query_id='q',
                    mapping_status='UNIQUE_EXACT_PLACEMENT',mapped_sequence_id='I',mapped_chromosome='1',
                    mapped_snp_bp=21,tag_start_bp=10,tag_end_bp=40,mapped_strand='+',assembly_base='A',
                    other_source_base='C',primary_chromosome=True)
    def marker(self,resource='Nacional_G7',key='k1'):
        return dict(resource=resource,source_row=2,marker_id=key,canonical_id=key,chromosome='1',position_bp=10)
    def test_linked_marker_explicit_scope(self):
        r,s=project_markers([self.marker()],[self.mapping()],'Diversity',['Nacional_G7'])
        self.assertTrue(r[0]['accepted_for_coordinate_comparison']);self.assertEqual(r[0]['coordinate_status'],'PROJECTED_VIA_UNIQUE_TAG_KEY')
        self.assertEqual(r[0]['genotype_allele_polarity'],'NOT_ESTABLISHED')
    def test_direct_source_row_not_wrong_name_join(self):
        old=self.marker('Diversity','wrong_name');r,_=project_markers([old],[self.mapping()],'Diversity',[])
        self.assertEqual(r[0]['reference_snp_bp'],21);self.assertEqual(r[0]['original_position_bp'],10)
    def test_direct_row_mismatch_fails(self):
        with self.assertRaises(DataError):project_markers([self.marker('Diversity','k1')],[self.mapping()],'Diversity',[])
    def test_amazonia_not_inferred_by_same_position(self):
        r,_=project_markers([self.marker('Amazonia')],[self.mapping()],'Diversity',['Nacional_G7'])
        self.assertFalse(r[0]['accepted_for_coordinate_comparison']);self.assertEqual(r[0]['coordinate_status'],'NO_SEQUENCE_IDENTITY_BRIDGE')
    def test_unknown_key_not_inferred(self):
        r,_=project_markers([self.marker(key='unknown')],[self.mapping()],'Diversity',['Nacional_G7'])
        self.assertEqual(r[0]['coordinate_status'],'NO_REFERENCE_TAG_KEY')
    def test_duplicate_reference_key_excluded(self):
        m=self.mapping();m2=dict(m,source_row=3)
        r,_=project_markers([self.marker()],[m,m2],'Diversity',['Nacional_G7'])
        self.assertEqual(r[0]['coordinate_status'],'AMBIGUOUS_REFERENCE_TAG_KEY')
    def test_multimap_excluded(self):
        m=dict(self.mapping(),mapping_status='MULTIPLE_EXACT_FULL_TAG_PLACEMENTS')
        r,_=project_markers([self.marker()],[m],'Diversity',['Nacional_G7'])
        self.assertFalse(r[0]['accepted_for_coordinate_comparison'])
    def test_duplicate_panel_key_excluded(self):
        m=self.marker();m2=dict(m,source_row=3)
        r,_=project_markers([m,m2],[self.mapping()],'Diversity',['Nacional_G7'])
        self.assertTrue(all(not a['accepted_for_coordinate_comparison'] for a in r))
    def test_unplaced_not_primary(self):
        m=dict(self.mapping(),primary_chromosome=False,mapped_chromosome=None,mapped_sequence_id='s1')
        r,_=project_markers([self.marker()],[m],'Diversity',['Nacional_G7'])
        self.assertFalse(r[0]['accepted_for_coordinate_comparison']);self.assertEqual(r[0]['reference_sequence_id'],'s1')
    def test_wrong_chromosome_not_fixed_in_original(self):
        old=dict(self.marker(),chromosome='7');r,_=project_markers([old],[self.mapping()],'Diversity',['Nacional_G7'])
        self.assertEqual(old['chromosome'],'7');self.assertEqual(r[0]['reference_chromosome'],'1');self.assertFalse(r[0]['source_chromosome_agrees'])
    def gene(self):
        return dict(gene_id='g1',gene_name='symbol1',chromosome='1',sequence_id='I',start_bp=10,end_bp=20,strand='+',biotype='protein_coding')
    def cand(self):return dict(gene_id='g1',chromosome='1',start_bp=10,end_bp=20,record_count=2,pathways='Sugar')
    def test_candidate_exact_interval(self):
        r=candidate_checks([self.cand()],[self.gene()],{})[0];self.assertTrue(r['accepted_original_interval'])
    def test_candidate_wrong_bounds_not_replaced(self):
        c=dict(self.cand(),start_bp=9);r=candidate_checks([c],[self.gene()],{})[0]
        self.assertFalse(r['accepted_original_interval']);self.assertEqual(c['start_bp'],9)
    def test_candidate_id_missing_not_nearest_gene(self):
        c=dict(self.cand(),gene_id='g2');r=candidate_checks([c],[self.gene()],{})[0]
        self.assertEqual(r['status'],'SOURCE_GENE_ID_NOT_FOUND_IN_THIS_ANNOTATION')
    def test_candidate_name_fallback_labeled(self):
        c=dict(self.cand(),gene_id='symbol1');r=candidate_checks([c],[self.gene()],{})[0]
        self.assertEqual(r['identity_field'],'GFF_Name_exact')
    def test_candidate_duplicate_ambiguous(self):
        r=candidate_checks([self.cand()],[self.gene(),self.gene()],{})[0]
        self.assertEqual(r['status'],'AMBIGUOUS_ANNOTATION_ID')
    def test_peak_context_does_not_rewrite_peak(self):
        old=self.marker();mapped,_=project_markers([old],[self.mapping()],'Diversity',['Nacional_G7'])
        p=dict(query_id='p1',chromosome='1',start_bp=10)
        r=peak_links([p],[old],mapped)[0]
        self.assertFalse(r['coordinate_rewrite_applied']);self.assertEqual(p['start_bp'],10)
        self.assertEqual(r['association_marker_identity'],'NOT_ESTABLISHED')
