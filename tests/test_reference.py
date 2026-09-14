import gzip
import io
import json
import logging
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from cacao_marker_harmonization.reference.common import DataError,sha256
from cacao_marker_harmonization.reference.reference import (iter_fasta, validate_sequence_header,validate_fasta_inventory,
                                      parse_attributes,read_gff,acquire_asset,download_asset)

class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.spec=dict(assembly_name='TEST_ASSEMBLY',assembly_accession='GCA_TEST',expected_total_bases=10,
                       chromosomes=['1'],chromosome_aliases={'I':'1'},gff_gene_features=['gene','ncRNA_gene','pseudogene'])
        self.log=logging.getLogger('test_ref')
    def tearDown(self):self.temp.cleanup()
    def fasta(self,s):
        p=self.root/'test.fa.gz'
        with gzip.open(p,'wt') as f:f.write(s)
        return p
    def test_stream_fasta_lowercase(self):
        r=list(iter_fasta(self.fasta('>I TEST_ASSEMBLY\nacgtNN\nACGT\n'),50))
        self.assertEqual(r[0][2],b'ACGTNNACGT')
    def test_duplicate_sequence_id(self):
        with self.assertRaises(DataError):list(iter_fasta(self.fasta('>I TEST_ASSEMBLY\nACGT\n>I TEST_ASSEMBLY\nACGT\n'),50))
    def test_empty_final_contig(self):
        with self.assertRaises(DataError):list(iter_fasta(self.fasta('>I TEST_ASSEMBLY\nACGT\n>II TEST_ASSEMBLY\n'),50))
    def test_missing_header(self):
        with self.assertRaises(DataError):list(iter_fasta(self.fasta('ACGT\n'),50))
    def test_invalid_sequence_char(self):
        with self.assertRaises(DataError):list(iter_fasta(self.fasta('>I TEST_ASSEMBLY\nACG!\n'),50))
    def test_contig_limit(self):
        with self.assertRaises(DataError):list(iter_fasta(self.fasta('>I TEST_ASSEMBLY\nACGTACGT\n'),7))
    def test_empty_fasta(self):
        with self.assertRaises(DataError):list(iter_fasta(self.fasta(''),50))
    def test_truncated_gzip(self):
        p=self.fasta('>I TEST_ASSEMBLY\nACGTACGT\n');p.write_bytes(p.read_bytes()[:-5])
        with self.assertRaises((EOFError,OSError)):list(iter_fasta(p,50))
    def test_assembly_header_required(self):
        with self.assertRaises(DataError):validate_sequence_header('I','I OTHER_ASSEMBLY',self.spec)
    def test_assembly_header_pass(self):validate_sequence_header('I','I chromosome:TEST_ASSEMBLY:I:1:10:1',self.spec)
    def test_fasta_total_length(self):
        self.assertEqual(validate_fasta_inventory([dict(chromosome='1',length_bp=8),dict(chromosome=None,length_bp=2)],self.spec),10)
    def test_partial_fasta_rejected(self):
        with self.assertRaises(DataError):validate_fasta_inventory([dict(chromosome='1',length_bp=8)],self.spec)
    def test_duplicate_primary_alias_rejected(self):
        with self.assertRaises(DataError):validate_fasta_inventory([dict(chromosome='1',length_bp=5),dict(chromosome='1',length_bp=5)],self.spec)
    def test_absent_primary_rejected(self):
        with self.assertRaises(DataError):validate_fasta_inventory([dict(chromosome=None,length_bp=10)],self.spec)
    def test_parse_percent_encoded_attributes(self):
        self.assertEqual(parse_attributes('ID=gene:g1;description=a%3Bb%20c;Name=g1')['description'],'a;b c')
    def test_gff_gene_features_only(self):
        p=self.root/'genes.gff3';p.write_text('##gff-version 3\n#!genome-build TEST_ASSEMBLY\nI\tx\tgene\t1\t8\t.\t+\t.\tID=gene:g1;Name=G1;biotype=protein_coding\nI\tx\tmRNA\t1\t8\t.\t+\t.\tID=transcript:g1;Parent=gene:g1\n', encoding='utf-8')
        gs,h,err,build=read_gff(p,{'I':10},{'I':'1'},self.spec)
        self.assertEqual(len(gs),1);self.assertTrue(build);self.assertEqual(gs[0]['gene_id'],'g1');self.assertFalse(err)
    def test_gff_out_of_bounds_not_repaired(self):
        p=self.root/'genes.gff3';p.write_text('##gff-version 3\nI\tx\tgene\t1\t8\t.\t+\t.\tID=g1\nI\tx\tgene\t8\t20\t.\t+\t.\tID=g2\n', encoding='utf-8')
        gs,h,err,build=read_gff(p,{'I':10},{'I':'1'},self.spec)
        self.assertEqual(len(gs),1);self.assertEqual(err[0]['issue'],'gene_outside_reference_sequence');self.assertFalse(build)
    def test_gff_duplicate_id_retained_for_ambiguity(self):
        p=self.root/'genes.gff3';p.write_text('##gff-version 3\nI\tx\tgene\t1\t5\t.\t+\t.\tID=g1\nI\tx\tgene\t6\t8\t.\t+\t.\tID=g1\n', encoding='utf-8')
        gs,_,err,_=read_gff(p,{'I':10},{'I':'1'},self.spec)
        self.assertEqual(len(gs),2);self.assertEqual(err[0]['issue'],'duplicate_gene_id')
    def test_gff_no_genes_fails(self):
        p=self.root/'genes.gff3';p.write_text('##gff-version 3\n', encoding='utf-8')
        with self.assertRaises(DataError):read_gff(p,{'I':10},{'I':'1'},self.spec)
    def spec_asset(self,p):
        return dict(self.spec,minimum_cache_free_gib=0,fasta=dict(filename=p.name,sha256=None,maximum_bytes=10000,urls=['https://example.org/ref.fa.gz']))
    def test_local_reference_no_download_no_copy(self):
        p=self.fasta('>I TEST_ASSEMBLY\nACGTACGTNN\n');cache=self.root/'cache'
        with patch('urllib.request.urlopen',side_effect=AssertionError('network not allowed')):
            a=acquire_asset('fasta',self.spec_asset(p),cache,[self.root],None,True,self.log)
        self.assertEqual(a.path,p);self.assertFalse(cache.exists())
    def test_offline_missing_fails(self):
        p=Path('test.fa.gz')
        with self.assertRaises(DataError):acquire_asset('fasta',self.spec_asset(p),self.root/'cache',[],None,True,self.log)
    def test_explicit_missing_fails(self):
        p=self.root/'missing.fa.gz'
        with self.assertRaises(DataError):acquire_asset('fasta',self.spec_asset(p),self.root/'cache',[],p,False,self.log)
    def test_cached_digest_change_fails(self):
        p=self.fasta('>I TEST_ASSEMBLY\nACGTACGTNN\n')
        (self.root/(p.name+'.provenance.json')).write_text(json.dumps(dict(sha256='bad')), encoding='utf-8')
        with self.assertRaises(DataError):acquire_asset('fasta',self.spec_asset(p),self.root,[],None,True,self.log)
    def test_configured_digest_change_fails(self):
        p=self.fasta('>I TEST_ASSEMBLY\nACGTACGTNN\n');s=self.spec_asset(p);s['fasta']['sha256']='bad'
        with self.assertRaises(DataError):acquire_asset('fasta',s,self.root,[],p,True,self.log)
    def response(self,data,url='https://example.org/ref.fa.gz',length=None):
        class Response(io.BytesIO):
            def __init__(self):
                super().__init__(data);self.headers={'Content-Length':str(len(data) if length is None else length)}
            def geturl(self):return url
        return Response()
    def test_http_download_mock_and_atomic_cache(self):
        data=gzip.compress(b'>I TEST_ASSEMBLY\nACGTACGTNN\n');dest=self.root/'download.fa.gz'
        with patch('urllib.request.urlopen',return_value=self.response(data)):
            url,digest=download_asset(['https://example.org/ref.fa.gz'],dest,'fasta',10000,self.log,attempts=1)
        self.assertEqual(sha256(dest),digest);self.assertEqual(dest.read_bytes(),data)
        self.assertFalse(list(self.root.glob('.download_*')))
    def test_html_download_rejected(self):
        dest=self.root/'download.fa.gz';data=gzip.compress(b'<html>Error</html>')
        with patch('urllib.request.urlopen',return_value=self.response(data)):
            with self.assertRaises(DataError):download_asset(['https://example.org/ref.fa.gz'],dest,'fasta',10000,self.log,attempts=1)
        self.assertFalse(dest.exists());self.assertFalse(list(self.root.glob('.download_*')))
    def test_incomplete_download_rejected(self):
        dest=self.root/'download.fa.gz';data=gzip.compress(b'>I TEST_ASSEMBLY\nACGT\n')
        with patch('urllib.request.urlopen',return_value=self.response(data,length=len(data)+2)):
            with self.assertRaises(DataError):download_asset(['https://example.org/ref.fa.gz'],dest,'fasta',10000,self.log,attempts=1)
        self.assertFalse(dest.exists())
    def test_non_https_download_rejected(self):
        with self.assertRaises(DataError):download_asset(['http://example.org/ref.fa.gz'],self.root/'d.fa.gz','fasta',10000,self.log,attempts=1)
    def test_download_size_bound(self):
        data=gzip.compress(b'>I TEST_ASSEMBLY\nACGT\n')
        with patch('urllib.request.urlopen',return_value=self.response(data)):
            with self.assertRaises(DataError):download_asset(['https://example.org/ref.fa.gz'],self.root/'d.fa.gz','fasta',5,self.log,attempts=1)
    def test_no_overwrite_reference(self):
        p=self.fasta('>I TEST_ASSEMBLY\nACGT\n');before=p.read_bytes()
        with patch('urllib.request.urlopen',return_value=self.response(before)):
            with self.assertRaises(DataError):download_asset(['https://example.org/ref.fa.gz'],p,'fasta',10000,self.log,attempts=1)
        self.assertEqual(p.read_bytes(),before)
