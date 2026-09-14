import copy,hashlib,json,logging,tempfile,unittest,zipfile
from pathlib import Path
import numpy as np
from cacao_marker_harmonization.coverage.common import DataError,sha256
from cacao_marker_harmonization.alignment.tags import source_tags,write_fasta,Tag
from cacao_marker_harmonization.alignment.refine import add_coding_query,project,background_sets,compare_previous,mapping_transitions
from cacao_marker_harmonization.coverage.numerics import Queries
from cacao_marker_harmonization.alignment.toolchain import reference_cache,Reference,validate_cache,safe_tool_extract
from cacao_marker_harmonization.configuration import load_study as load_config

class ReferenceCacheTests(unittest.TestCase):
    def setUp(self):self.log=logging.getLogger('test')
    def fixture(self,d):
        p=Path(d)/'genome.fa';seq='ACTGNNACGT';p.write_text('>I assembly\n'+seq+'\n', encoding='utf-8')
        return p,{'sha256':sha256(p)},[dict(sequence_id='I',chromosome='1',length_bp=len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())]
    def test_readonly_reference(self):
        with tempfile.TemporaryDirectory() as d:
            p,a,i=self.fixture(d);cache=Path(d)/'cache';h=sha256(p);f,j=reference_cache(p,a,i,cache,self.log)
            with Reference(f,j) as ref:self.assertEqual(ref.slice('I',2,5),b'TGN')
            self.assertEqual(h,sha256(p));self.assertEqual(reference_cache(p,a,i,cache,self.log),(f,j))
    def test_inventory_hash_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            p,a,i=self.fixture(d);i[0]['sequence_sha256']='0'*64
            with self.assertRaises(DataError):reference_cache(p,a,i,Path(d)/'cache',self.log)
    def test_cache_tampering(self):
        with tempfile.TemporaryDirectory() as d:
            p,a,i=self.fixture(d);cache=Path(d)/'cache';f,j=reference_cache(p,a,i,cache,self.log);f.write_text('bad', encoding='utf-8')
            with self.assertRaises(DataError):reference_cache(p,a,i,cache,self.log)
    def test_slice_bounds(self):
        with tempfile.TemporaryDirectory() as d:
            p,a,i=self.fixture(d);f,j=reference_cache(p,a,i,Path(d)/'cache',self.log)
            with Reference(f,j) as ref:
                for args in [('I',-1,2),('I',0,100),('missing',0,1)]:
                    with self.assertRaises(DataError):ref.slice(*args)
    def test_unsafe_tool_zip(self):
        with tempfile.TemporaryDirectory() as d:
            z=Path(d)/'x.zip';dest=Path(d)/'out';dest.mkdir()
            with zipfile.ZipFile(z,'w') as f:f.writestr('../escape','x')
            with self.assertRaises(DataError):safe_tool_extract(z,dest)
    def test_tool_symlink_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            z=Path(d)/'x.zip';dest=Path(d)/'out';dest.mkdir();m=zipfile.ZipInfo('link');m.external_attr=(0o120777<<16)
            with zipfile.ZipFile(z,'w') as f:f.writestr(m,'/etc')
            with self.assertRaises(DataError):safe_tool_extract(z,dest)
    def test_tool_duplicate_path_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            z=Path(d)/'x.zip';dest=Path(d)/'out';dest.mkdir()
            with zipfile.ZipFile(z,'w') as f:f.writestr('A','1');f.writestr('a','2')
            with self.assertRaises(DataError):safe_tool_extract(z,dest)

class TagTests(unittest.TestCase):
    def data(self):
        seq='A'*30+'C'+'G'*30
        e=dict(source_row=2,marker_id='mk',reference_key='key',five_flank='A'*30,three_flank='G'*30,variation='C>T',reference_sequence=seq,ref_lookup_status='exact_ref_tag')
        p=dict(query_id='diversity_row_2',source_marker_id='mk',reference_key='key',eligibility='ELIGIBLE',mapping_status='NO_EXACT_FULL_TAG_MATCH',sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),source_snp_index_0based=30,tag_length=61)
        return e,p
    def test_both_tags(self):
        e,p=self.data();t,a=source_tags([e],[p]);self.assertEqual(t['diversity_row_2'].alt[30],'T')
    def test_evidence_change(self):
        e,p=self.data();e['five_flank']='T'+e['five_flank'][1:]
        with self.assertRaises(DataError):source_tags([e],[p])
    def test_offset_change(self):
        e,p=self.data();p['source_snp_index_0based']=29
        with self.assertRaises(DataError):source_tags([e],[p])
    def test_all_prior_tags_required(self):
        e,p=self.data()
        with self.assertRaises(DataError):source_tags([], [p])
    def test_duplicate_source_row(self):
        e,p=self.data()
        with self.assertRaises(DataError):source_tags([e,e],[p])
    def test_fasta_no_overwrite(self):
        e,p=self.data();t,_=source_tags([e],[p])
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/'tags.fa';write_fasta(f,t)
            self.assertEqual(f.read_text(encoding='utf-8').count('>'),2)
            with self.assertRaises(FileExistsError):write_fasta(f,t)

class CodingTests(unittest.TestCase):
    def setUp(self):
        self.genes=[dict(gene_id='g'+str(i),query_id='g'+str(i),chromosome='1',start_bp=100*i+1,end_bp=100*i+30,biotype='protein_coding') for i in range(12)]
        self.genes[1]['biotype']='nontranslating_CDS'
        self.q={'unique_candidate_genes':Queries('unique_candidate_genes',self.genes[:2])}
    def test_full_set_preserved(self):
        q,l=add_coding_query(self.q);self.assertEqual(len(q['unique_candidate_genes']),2);self.assertEqual(len(q['protein_coding_candidate_genes']),1);self.assertEqual(len(l),2)
    def test_empty_coding_blocks(self):
        self.genes[0]['biotype']='nontranslating_CDS'
        with self.assertRaises(DataError):add_coding_query(self.q)
    def test_gene_density_deficit_still_blocks(self):
        # Give the rare biotype a noncandidate donor, but then deliberately exclude
        # all donors of coding type; the entire model is blocked, not reweighted.
        self.genes[2]['biotype']='nontranslating_CDS';q,_=add_coding_query(self.q)
        cfg={'seed':1729,'comparison_resources':['r'],'windows_bp':[0,50],
             'background':{'enabled':True,'replicates':5,'length_bins':1,'local_gene_density_bins':1,'local_gene_density_radius_bp':100,
             'models':['chromosome_biotype_length'],'query_set':'unique_candidate_genes','query_sets':['unique_candidate_genes','protein_coding_candidate_genes'],'formal_p_values':False}}
        ps={'r':{arm:{'1':np.array([1,101,301])} for arm in ['A','B','C']}}
        out=background_sets(self.genes,{'g0','g1'},q,{'1':2000},ps,cfg,logging.getLogger('t'))
        self.assertEqual({s['n_target_genes'] for s in out['status']},{1,2})
        self.assertTrue(all(s['status']=='COMPLETED_DESCRIPTIVE_BACKGROUND' for s in out['status']))
        self.assertTrue(all(r['excluded_candidate'] for r in out['background_gene_distances'] if r['gene_id'] in ('g0','g1')))

class ProjectionTests(unittest.TestCase):
    def setUp(self):
        self.cfg={'comparison_resources':['Diversity','Nacional_G7','Nacional_MAF5']}
        self.old=[dict(resource='Diversity',source_row=2,marker_id='mk',canonical_id='mk',chromosome='1',position_bp=50),dict(resource='Nacional_G7',source_row=2,marker_id='key',canonical_id='key',chromosome='1',position_bp=20)]
        self.m=dict(query_id='diversity_row_2',source_row=2,source_marker_id='mk',reference_key='key',placement_status='RETAINED_PRIMARY_CHROMOSOME',accepted_for_mapping=True,mapped_sequence_id='I',mapped_chromosome='1',mapped_snp_bp=50,tag_start_bp=20,tag_end_bp=88,mapped_strand='+',assembly_base='A',other_source_base='C',ref_cigar='69M',alt_cigar='69M')
    def test_identity_projection_labeled(self):
        r=project(self.old,[self.m],self.cfg);self.assertTrue(r[1]['accepted_for_coordinate_comparison']);self.assertIn('not_independent',r[1]['identity_basis'])
    def test_ambiguous_source_key(self):
        x=dict(self.m,query_id='diversity_row_3',source_row=3);r=project(self.old,[self.m,x],self.cfg);self.assertFalse(r[1]['accepted_for_coordinate_comparison'])
    def test_context_not_projected(self):
        o=dict(self.old[1],resource='Amazonia');r=project([o],[self.m],self.cfg);self.assertFalse(r[0]['accepted_for_coordinate_comparison'])
    def test_same_marker_key_twice(self):
        old=self.old+[dict(self.old[1],source_row=3)];r=project(old,[self.m],self.cfg);self.assertFalse(r[1]['accepted_for_coordinate_comparison'])
    def test_no_automatic_prior_acceptance(self):
        prev=dict(query_id='diversity_row_2',mapping_status='UNIQUE_EXACT_PLACEMENT',primary_chromosome=True,mapped_sequence_id='I',mapped_snp_bp=50,mapped_strand='+')
        cur=dict(self.m,accepted_for_mapping=False,placement_status='INSUFFICIENT_SCORE_SEPARATION');r=mapping_transitions([prev],[cur]);self.assertEqual(r[0]['transition'],'EXACT_ONLY_RETAINED_NOW_EXCLUDED')
    def test_original_arm_drift_error(self):
        a=dict(resource='r',arm='A',query_set='g',flank_bp=0,n_queries=2,n_recovered=1,n_marker_records=2,fraction_recovered=.5)
        with self.assertRaises(DataError):compare_previous([a],[dict(a,n_recovered=2)])
    def test_query_denominator_drift_error(self):
        a=dict(resource='r',arm='B',query_set='g',flank_bp=0,n_queries=2,n_recovered=1,n_marker_records=2,fraction_recovered=.5)
        with self.assertRaises(DataError):compare_previous([a],[dict(a,n_queries=3)])

class ConfigTests(unittest.TestCase):
    def test_complete_config(self):self.assertEqual(load_config(Path(__file__).resolve().parents[1]/'config/study.json')['alignment']['version'],'2.5.5')
    def test_bounded_threads(self):
        cfg=json.loads((Path(__file__).resolve().parents[1]/'config/study.json').read_text(encoding='utf-8'));cfg['alignment']['threads']=100
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'c.json';p.write_text(json.dumps(cfg), encoding='utf-8')
            with self.assertRaises(DataError):load_config(p)
