"""Boundary checks for source genotype interpretation and missingness."""
import unittest
import numpy as np
from genotype_audit import Matrix, Panel, canonical_marker, decode_matrix, sample_symbol_qc, compare_panel_cells, marker_symbol_qc

class GenotypeAuditTests(unittest.TestCase):
    def test_canonical_ids_keep_snp_offset(self):
        self.assertEqual(canonical_marker('123.F.0.7'),'123|F|0--7')
        self.assertNotEqual(canonical_marker('123.F.0.7'),canonical_marker('123.F.0.8'))

    def test_missing_and_unsupported_calls_not_zero_filled(self):
        m=Matrix('x',['a','b','c','d'],[{'source_row':2,'marker_id':'m','alleles_reported':'A/C'}],[('A','M','C','N')],'s','test',[])
        cfg={'encoding':{'symbol_to_dosage':{'A':0,'M':1,'C':2},'missing_symbols':['N'],'declared_allele_set':['A','C']}}
        d,state,_,bad=decode_matrix(m,cfg)
        np.testing.assert_equal(d[0,:3],[0,1,2]);self.assertTrue(np.isnan(d[0,3]));self.assertEqual(len(bad),0)
        m.calls=[('A','X','C','N')]
        d,state,_,bad=decode_matrix(m,cfg)
        self.assertTrue(np.isnan(d).all());self.assertEqual(state.tolist(),['UNSUPPORTED_SYMBOLS']);self.assertEqual(len(bad),1)

    def test_other_allele_declaration_not_interpreted_as_acm(self):
        m=Matrix('x',['a'],[{'source_row':2,'marker_id':'m','alleles_reported':'G/T'}],[('M',)],'s','test',[])
        cfg={'encoding':{'symbol_to_dosage':{'A':0,'M':1,'C':2},'missing_symbols':['N'],'declared_allele_set':['A','C']}}
        d,state,_,bad=decode_matrix(m,cfg)
        self.assertTrue(np.isnan(d).all());self.assertEqual(state.tolist(),['UNSUPPORTED_DECLARED_ALLELES'])

    def test_literal_m_fraction_uses_called_denominator(self):
        m=Matrix('x',['a'],[{},{}],[('M',),('N',)],'s','test',[])
        row=sample_symbol_qc(m,[0,1],'t',{'N'},set())[0]
        self.assertEqual(row['call_rate'],0.5);self.assertEqual(row['M_fraction_of_called'],1.0)
        self.assertIn('Conditional',row['interpretation'])

    def test_missingness_discrepancy_is_not_equal(self):
        markers=[{'marker_key':'m1'},{'marker_key':'m2'}]
        a=Panel('a',['s'],markers,np.array([[0.],[np.nan]]),{},{});
        b=Panel('b',['s'],markers,np.array([[np.nan],[np.nan]]),{},{});
        row=compare_panel_cells(a,b,'t')[0]
        self.assertEqual(row['one_missing_calls'],1);self.assertEqual(row['both_missing_calls'],1)
        self.assertFalse(row['exact_equal_including_missing'])

    def test_conditional_marker_frequencies(self):
        marker={'source_row':2,'marker_id':'m','canonical_id':'m','alleles_reported':'A/C'}
        m=Matrix('x',['s1','s2','s3','s4'],[marker],[('A','M','C','N')],'s','test',[])
        row=marker_symbol_qc(m,{'N'})[0]
        self.assertAlmostEqual(row['conditional_heterozygous_fraction'],1/3)
        self.assertEqual(row['conditional_C_state_allele_frequency'],0.5)
        self.assertEqual(row['conditional_two_p_one_minus_p'],0.5)
        self.assertFalse(row['HWE_test_performed'])

if __name__=='__main__':unittest.main()
