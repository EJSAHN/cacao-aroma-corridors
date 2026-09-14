import unittest,tempfile,json,logging
from pathlib import Path
import numpy as np
from cacao_marker_harmonization.io.export import write_workbook,cell_xml,clean_value,column_name
from cacao_marker_harmonization.coverage.inputs import read_sheet
from cacao_marker_harmonization.coverage.common import DataError
from cacao_marker_harmonization.coverage.numerics import Queries
from cacao_marker_harmonization.coverage.analysis import background_diagnostics

class WorkbookTests(unittest.TestCase):
    def test_column_names(self):self.assertEqual([column_name(n) for n in [1,26,27,703]],['A','Z','AA','AAA'])
    def test_roundtrip(self):
        rows=[dict(text='a<&"é 中文',i=123,f=.867761260261815,flag=True,empty=None),dict(text='=1+1',i=-1,f=0.,flag=False,empty='NA')]
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.xlsx';write_workbook(p,{'data':rows},[('note','values only')]);self.assertEqual(read_sheet(p,'data'),rows)
    def test_empty_sheet(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.xlsx';write_workbook(p,{'data':[]},[]);self.assertEqual(read_sheet(p,'data'),[])
    def test_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.xlsx';write_workbook(p,{'data':[]},[])
            with self.assertRaises(FileExistsError):write_workbook(p,{'data':[]},[])
    def test_formula_plain_text(self):self.assertNotIn('<f>',cell_xml('A1','=SUM(B1:B2)'))
    def test_nonfinite_rejected(self):
        for v in [float('nan'),float('inf')]:
            with self.assertRaises(ValueError):clean_value(v)
    def test_text_not_truncated(self):
        with self.assertRaises(ValueError):clean_value('x'*32768)
    def test_bad_sheet_name(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):write_workbook(Path(d)/'x.xlsx',{'bad/name':[]},[])
    def test_case_duplicate_sheet(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):write_workbook(Path(d)/'x.xlsx',{'Data':[],'data':[]},[])

class BackgroundTests(unittest.TestCase):
    def setUp(self):
        self.genes=[dict(gene_id='g'+str(i),chromosome='1',start_bp=100*i+1,end_bp=100*i+10,biotype='coding') for i in range(12)]
        self.q=Queries('unique_candidate_genes',[dict(self.genes[0],query_id='g0'),dict(self.genes[1],query_id='g1')])
        self.cfg={'seed':1729,'comparison_resources':['r'],'windows_bp':[0,10,50],
            'background':{'enabled':True,'replicates':10,'length_bins':1,'local_gene_density_bins':1,'local_gene_density_radius_bp':250,
            'models':['chromosome_biotype_length'],'query_set':'unique_candidate_genes','formal_p_values':False}}
        self.ps={'r':{arm:{'1':np.array([2,202,602,902])} for arm in ('A','B','C')}}
    def test_completes_without_p_values(self):
        res=background_diagnostics(self.genes,{'g0','g1'},{'unique_candidate_genes':self.q},{'1':2000},self.ps,self.cfg,logging.getLogger('test'))
        self.assertEqual(res['status'][0]['status'],'COMPLETED_DESCRIPTIVE_BACKGROUND')
        self.assertFalse(any('p_value' in k for r in res['summary'] for k in r))
    def test_same_draw_all_arms(self):
        res=background_diagnostics(self.genes,{'g0','g1'},{'unique_candidate_genes':self.q},{'1':2000},self.ps,self.cfg,logging.getLogger('test'))
        rows=res['replicate_counts'];A=[r['n_recovered_0_bp'] for r in rows if r['arm']=='A'];C=[r['n_recovered_0_bp'] for r in rows if r['arm']=='C'];self.assertEqual(A,C)
    def test_excludes_all_candidates_from_pool(self):
        res=background_diagnostics(self.genes,{'g0','g1','g2'},{'unique_candidate_genes':self.q},{'1':2000},self.ps,self.cfg,logging.getLogger('test'))
        self.assertEqual(res['stratum_requirements'][0]['available_noncandidate_genes'],9)
    def test_deficit_not_partial_target(self):
        # One background gene available for two target genes: block the whole model.
        res=background_diagnostics(self.genes,{'g'+str(i) for i in range(11)},{'unique_candidate_genes':self.q},{'1':2000},self.ps,self.cfg,logging.getLogger('test'))
        self.assertEqual(res['status'][0]['status'],'BLOCKED_INSUFFICIENT_STRATA');self.assertEqual(res['summary'],[])
    def test_gene_id_duplicate_error(self):
        with self.assertRaises(DataError):background_diagnostics(self.genes+[self.genes[0]],{'g0','g1'},{'unique_candidate_genes':self.q},{'1':2000},self.ps,self.cfg,logging.getLogger('test'))
    def test_interval_mismatch_error(self):
        q=Queries('unique_candidate_genes',[dict(self.genes[0],query_id='g0',end_bp=12)])
        with self.assertRaises(DataError):background_diagnostics(self.genes,{'g0','g1'},{'unique_candidate_genes':q},{'1':2000},self.ps,self.cfg,logging.getLogger('test'))
    def test_background_deterministic(self):
        a=background_diagnostics(self.genes,{'g0','g1'},{'unique_candidate_genes':self.q},{'1':2000},self.ps,self.cfg,logging.getLogger('test'))
        b=background_diagnostics(self.genes,{'g0','g1'},{'unique_candidate_genes':self.q},{'1':2000},self.ps,self.cfg,logging.getLogger('test'))
        self.assertEqual(a,b)
