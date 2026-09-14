import unittest
from cacao_marker_harmonization.coverage.inputs import build_queries
from cacao_marker_harmonization.coverage.common import DataError


def gene(g='g',start=10,end=20):return dict(gene_id=g,chromosome='1',start_bp=start,end_bp=end,pathways='P',original_rows='3')
def anno(g='g',start=10,end=20):return dict(source_gene_id=g,source_chromosome='1',source_start_bp=start,source_end_bp=end,
    accepted_original_interval=True,status='EXACT_ID_AND_INTERVAL',annotation_chromosome='1',annotation_start_bp=start,annotation_end_bp=end,biotype='coding')

class QueriesTests(unittest.TestCase):
    def test_preserves_candidate_ID_and_bounds(self):
        q,a,ids=build_queries([gene()],[anno()],[dict(gene_id='g',pathway='P')],{}, {'1':100})
        self.assertEqual(q['unique_candidate_genes'].rows[0]['start_bp'],10);self.assertEqual(ids,{'g'})
    def test_membership_duplicate_not_double_gene(self):
        q,_,_=build_queries([gene()],[anno()],[dict(gene_id='g',pathway='P')]*2,{}, {'1':100})
        self.assertEqual(len(q['pathway:P']),1)
    def test_undefined_membership_rejected(self):
        with self.assertRaises(DataError):build_queries([gene()],[anno()],[dict(gene_id='missing',pathway='P')],{}, {'1':100})
    def test_source_peak_not_shifted(self):
        source={'master_unique_peaks':[dict(query_id='p',chromosome='1',start_bp=15,end_bp=15)]}
        q,_,_=build_queries([gene()],[anno()],[],source, {'1':100})
        self.assertEqual(q['master_unique_peaks'].rows[0]['start_bp'],15)
    def test_source_interval_outside_bounds_excluded(self):
        source={'intervals':[dict(query_id='p',chromosome='1',start_bp=15,end_bp=115)]}
        q,a,_=build_queries([gene()],[anno()],[],source, {'1':100})
        self.assertNotIn('intervals',q);self.assertFalse(a[-1]['included'])
    def test_duplicate_peak_coordinates_error(self):
        source={'p':[dict(query_id=i,chromosome='1',start_bp=15,end_bp=15) for i in ['p1','p2']]}
        with self.assertRaises(DataError):build_queries([gene()],[anno()],[],source, {'1':100})
    def test_candidate_annotation_conflict_error(self):
        a=anno();a['annotation_end_bp']=21
        with self.assertRaises(DataError):build_queries([gene()],[a],[],{}, {'1':100})
    def test_candidate_manifest_mismatch_error(self):
        with self.assertRaises(DataError):build_queries([gene()],[anno('other')],[],{}, {'1':100})
    def test_source_metadata_changed_error(self):
        with self.assertRaises(DataError):build_queries([gene()],[anno(end=21)],[],{}, {'1':100})
