import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from cacao_marker_harmonization.ld.common import DataError,canonical_marker
from cacao_marker_harmonization.ld.readers import Matrix,read_matrix
from cacao_marker_harmonization.ld.genotypes import decode_matrix,row_qc,prepare_panel,prepare_common,raw_concordance
from cacao_marker_harmonization.io.export import write_workbook

def config():
    c=json.loads((Path(__file__).parents[1]/'config/study.json').read_text(encoding='utf-8'))['dosage_ld']
    c['sensitivity_sample_ids']=['s0'];c['qc']['min_pair_samples']=3;c['qc']['min_markers']=3
    return c

def fixture(name='x'):
    samples=['s0','s1','s2','s3','s4','s5']
    calls=[('C','A','M','C','M','A'),('A','C','M','A','M','C'),('M','A','C','C','A','M'),('A','M','C','M','A','C')]
    markers=[];ledger=[]
    for i in range(4):
        mid=f'{i+1}|F|0--5'
        markers.append(dict(resource=name,source_file='test.xls',source_sheet='data',source_row=i+2,marker_id=mid,canonical_id=mid,
                       chromosome='1',position_bp=i*100+1,alleles_reported='A/C',strand_reported='+',record_index=i))
        ledger.append(dict(resource=name,record_key=f'{name}:row:{i+2}',original_row=i+2,marker_id=mid,in_B=True,in_C=True,
                    original_chromosome='1',original_position_bp=i*100+1,reference_chromosome='1',reference_snp_bp=i*100+6,
                    source_tag_key=mid,identity_basis='key_link'))
    return Matrix(name,samples,markers,calls,'data','test',[]),ledger

class GenotypeTests(unittest.TestCase):
    def test_acm(self):
        m,_=fixture();d,s,_,_=decode_matrix(m,config());np.testing.assert_equal(d[0],[2,0,1,2,1,0]);self.assertEqual(set(s),{'OK'})
    def test_missing_is_nan(self):
        m,_=fixture();m.calls[0]=('N','A','M','C','M','A');d,_,_,_=decode_matrix(m,config());self.assertTrue(np.isnan(d[0,0]))
    def test_gtm_not_reinterpreted(self):
        m,_=fixture();m.calls[0]=('G','T','M','G','T','M');d,s,_,bad=decode_matrix(m,config());self.assertEqual(s[0],'UNSUPPORTED_SYMBOLS');self.assertTrue(np.isnan(d[0]).all());self.assertEqual(len(bad),4)
    def test_declared_alleles_must_match(self):
        m,_=fixture();m.markers[0]['alleles_reported']='G/T';_,s,_,_=decode_matrix(m,config());self.assertEqual(s[0],'UNSUPPORTED_DECLARED_ALLELES')
    def test_maf_pair_convention(self):
        c=config();d=np.array([[0,1,2,1.],[0,0,0,0],[0,2,np.nan,np.nan]])
        q=row_qc(d,c['qc']);self.assertEqual(q['maf'][0],.5);self.assertFalse(q['keep'][1]);self.assertFalse(q['keep'][2])
    def test_fixed_qc_intersection(self):
        c=config();m,l=fixture();m.calls[0]=('C','A','A','A','A','A')
        p,r=prepare_panel(m,l,c);self.assertEqual(len(p.markers),3);self.assertFalse(r['marker_qc'][0]['retained_for_fixed_ld'])
    def test_selected_samples_stay_in_primary(self):
        m,l=fixture();p,_=prepare_panel(m,l,config());self.assertEqual(len(p.groups['all_qc_samples']),6);self.assertEqual(len(p.groups['omit_flagged_samples']),5)
    def test_no_call_edits(self):
        m,l=fixture();before=copy.deepcopy(m.calls);prepare_panel(m,l,config());self.assertEqual(before,m.calls)
    def test_missing_sample_is_error(self):
        m,l=fixture();c=config();c['sensitivity_sample_ids']=['absent']
        with self.assertRaises(DataError):prepare_panel(m,l,c)
    def test_source_marker_drift(self):
        m,l=fixture();l[0]['marker_id']='wrong'
        with self.assertRaises(DataError):prepare_panel(m,l,config())
    def test_source_position_drift(self):
        m,l=fixture();l[0]['original_position_bp']=99
        with self.assertRaises(DataError):prepare_panel(m,l,config())
    def test_source_row_drift(self):
        m,l=fixture();l.pop()
        with self.assertRaises(DataError):prepare_panel(m,l,config())
    def test_duplicate_canonical_is_error(self):
        m,l=fixture();m.markers[0]['canonical_id']=m.markers[1]['canonical_id']
        with self.assertRaises(DataError):prepare_panel(m,l,config())
    def test_duplicate_coordinate_excludes_all(self):
        m,l=fixture();c=config();l[0]['reference_snp_bp']=l[1]['reference_snp_bp'];c['qc']['min_markers']=2
        p,r=prepare_panel(m,l,c);self.assertEqual(len(p.markers),2);self.assertEqual(sum(x['duplicate_corrected_position'] for x in r['marker_qc']),2)
    def test_unplaced_is_not_promoted(self):
        m,l=fixture();l[0]['in_B']=False;l[0]['in_C']=False;l[0]['reference_chromosome']='scaffold_1';l[0]['reference_snp_bp']=None
        p,r=prepare_panel(m,l,config());self.assertEqual(len(p.markers),3)
    def test_accepted_nonprimary_is_error(self):
        m,l=fixture();l[0]['reference_chromosome']='scaffold_1'
        with self.assertRaises(DataError):prepare_panel(m,l,config())
    def test_duplicate_corrected_does_not_drop_original_collisions(self):
        m,l=fixture();m.markers[0]['position_bp']=101;l[0]['original_position_bp']=101
        p,_=prepare_panel(m,l,config());self.assertEqual(len(p.markers),4)
    def test_dart_name_alias_is_explicit(self):
        self.assertEqual(canonical_marker('123.F.0.2'),'123|F|0--2');self.assertNotEqual(canonical_marker('sample.a'),canonical_marker('sample_a'))
    def test_excel_error_fails(self):
        m,l=fixture();m.errors=[{'error':'#REF!'}]
        with self.assertRaises(DataError):prepare_panel(m,l,config())

class SharedTests(unittest.TestCase):
    def test_concordance_flags_not_auto_changed(self):
        a,_=fixture('a');b,_=fixture('b');b.calls[0]=('A','C','M','C','M','A')
        r=raw_concordance(a,b,config());self.assertEqual(r['summary'][0]['discordant_calls'],2);self.assertEqual(r['summary'][0]['discordance_outside_declared_samples'],1)
    def test_same_samples_are_aligned_by_id(self):
        a,_=fixture('a');b,_=fixture('b');b.samples=b.samples[::-1];b.calls=[x[::-1] for x in b.calls]
        self.assertEqual(raw_concordance(a,b,config())['summary'][0]['discordant_calls'],0)
    def test_marker_order_not_assumed(self):
        a,_=fixture('a');b,_=fixture('b');b.markers=b.markers[::-1];b.calls=b.calls[::-1]
        self.assertEqual(raw_concordance(a,b,config())['summary'][0]['discordant_calls'],0)
    def test_common_markers_fixed(self):
        c=config();a,la=fixture('a');b,lb=fixture('b');a,_=prepare_panel(a,la,c);b,_=prepare_panel(b,lb,c)
        a,b,audit,shared,omit=prepare_common(a,b,c)
        self.assertEqual(len(shared),6);self.assertEqual(len(omit),5);self.assertEqual(len(a.markers),4)
        self.assertEqual([x['marker_key'] for x in a.markers],[x['marker_key'] for x in b.markers])
    def test_common_coordinate_conflict_not_merged(self):
        c=config();a,la=fixture('a');b,lb=fixture('b');lb[0]['reference_snp_bp']=50
        a,_=prepare_panel(a,la,c);b,_=prepare_panel(b,lb,c)
        aa,bb,audit,_,_=prepare_common(a,b,c)
        self.assertEqual(len(aa.markers),3);self.assertEqual(audit[0]['reason'],'COORDINATE_OR_IDENTITY_CONFLICT')

class ReaderTests(unittest.TestCase):
    def test_marker_rows_not_sample_rows(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'matrix.xlsx';write_workbook(p,{'data':[{'rs#':'1|F|0--2','alleles':'A/C','chrom':1,'pos':3,'s0':'A','s1':'M','s2':'C'}]},[])
            m=read_matrix(p,'x',{'matrix_sheet':'data'});self.assertEqual(m.samples,['s0','s1','s2']);self.assertEqual(m.calls,[('A','M','C')])
    def test_multiple_sheets_need_selection(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'matrix.xlsx';write_workbook(p,{'data':[{'x':1}]},[])
            with self.assertRaises(DataError):read_matrix(p,'x',{})
