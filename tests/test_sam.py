import copy,gzip,json,tempfile,unittest
from pathlib import Path
from dataclasses import replace
from cacao_marker_harmonization.alignment.tags import Tag,revcomp
from cacao_marker_harmonization.alignment.sam import cigar_ops,parse_record,allele_decision,combine_alleles,parse_sam,Placement
from cacao_marker_harmonization.coverage.common import DataError

CFG=json.loads((Path(__file__).resolve().parents[1]/'config/study.json').read_text(encoding='utf-8'))['alignment']
class MemoryRef:
    def __init__(self,seqs):self.seqs=seqs;self.entries={k:{'length':len(v)} for k,v in seqs.items()}
    def slice(self,k,a,b):
        if k not in self.seqs or not 0<=a<=b<=len(self.seqs[k]):raise DataError('bounds')
        return self.seqs[k][a:b].encode()

def samline(tag,allele='R',pos=11,flag=0,cigar=None,nm=0,score=0,seqid='I',xs=None):
    seq=tag.seq(allele)
    if flag&16:seq=revcomp(seq)
    f=[tag.query_id+'__'+allele,str(flag),seqid,str(pos),'42',cigar or str(len(seq))+'M','*','0','0',seq,'*','AS:i:'+str(score),'NM:i:'+str(nm)]
    if xs is not None:f.append('XS:i:'+str(xs))
    return '\t'.join(f)+'\n'

def phit(score=0,pos=10,snp=20,reason='ELIGIBLE_PLACEMENT',xs=None):
    return Placement('q__R','I',pos,pos+49,'+','50M',score,0,42,snp,'A','C',0,0,None,reason,xs)

class SamNumerics(unittest.TestCase):
    def setUp(self):
        seq='ACGTTGCACTTAGCGACATCGATGTACGATTCGGAATCTAGCTACGTAG'
        self.t=Tag('q',2,'key','mk',seq,seq[:20]+'A'+seq[21:],20)
        assert self.t.ref[20]!='A'
    def ref(self,seq=None):return MemoryRef({'I':'N'*10+(seq or self.t.ref)+'N'*20})
    def test_forward(self):
        p=parse_record(samline(self.t),self.t,'R',self.ref(),CFG);self.assertEqual(p.snp,31);self.assertEqual(p.reason,'ELIGIBLE_PLACEMENT')
    def test_reverse(self):
        p=parse_record(samline(self.t,flag=16),self.t,'R',self.ref(revcomp(self.t.ref)),CFG);self.assertEqual(p.snp,11+len(self.t.ref)-1-20);self.assertEqual(p.strand,'-')
    def test_focal_mismatch_not_nonfocal(self):
        p=parse_record(samline(self.t,'A',nm=1,score=-6),self.t,'A',self.ref(),CFG);self.assertEqual(p.nonfocal_edits,0);self.assertEqual(p.reason,'ELIGIBLE_PLACEMENT')
    def test_reverse_alt(self):
        p=parse_record(samline(self.t,'A',flag=16,nm=1,score=-6),self.t,'A',self.ref(revcomp(self.t.ref)),CFG);self.assertEqual(p.nonfocal_edits,0)
    def test_source_sequence_mismatch(self):
        s=samline(self.t).replace(self.t.ref,'A'*len(self.t.ref))
        with self.assertRaises(DataError):parse_record(s,self.t,'R',self.ref(),CFG)
    def test_nm_mismatch(self):
        with self.assertRaises(DataError):parse_record(samline(self.t,nm=1),self.t,'R',self.ref(),CFG)
    def test_wrong_query_length(self):
        with self.assertRaises(DataError):parse_record(samline(self.t,cigar='5M'),self.t,'R',self.ref(),CFG)
    def test_unmapped(self):self.assertIsNone(parse_record(samline(self.t,flag=4),self.t,'R',self.ref(),CFG))
    def test_paired_rejected(self):
        with self.assertRaises(DataError):parse_record(samline(self.t,flag=1),self.t,'R',self.ref(),CFG)
    def test_supplementary_rejected(self):
        with self.assertRaises(DataError):parse_record(samline(self.t,flag=2048),self.t,'R',self.ref(),CFG)
    def test_secondary_255_not_used(self):
        s=samline(self.t,flag=256).replace('\t42\t','\t255\t');p=parse_record(s,self.t,'R',self.ref(),CFG);self.assertEqual(p.mapq,255)
    def test_insertion_distal(self):
        seq=self.t.ref[:5]+self.t.ref[6:];c='5M1I'+str(len(self.t.ref)-6)+'M'
        p=parse_record(samline(self.t,cigar=c,nm=1,score=-8),self.t,'R',self.ref(seq),CFG);self.assertEqual(p.snp,30);self.assertEqual(p.gap_bases,1);self.assertEqual(p.reason,'ELIGIBLE_PLACEMENT')
    def test_deletion_distal(self):
        seq=self.t.ref[:5]+'G'+self.t.ref[5:];c='5M1D'+str(len(self.t.ref)-5)+'M'
        p=parse_record(samline(self.t,cigar=c,nm=1,score=-8),self.t,'R',self.ref(seq),CFG);self.assertEqual(p.snp,32);self.assertEqual(p.reason,'ELIGIBLE_PLACEMENT')
    def test_focal_insertion_blocks(self):
        seq=self.t.ref[:20]+self.t.ref[21:];c='20M1I'+str(len(self.t.ref)-21)+'M'
        p=parse_record(samline(self.t,cigar=c,nm=1,score=-8),self.t,'R',self.ref(seq),CFG);self.assertEqual(p.reason,'FOCAL_SNP_NOT_ON_REFERENCE_BASE')
    def test_gap_near_snp_blocks(self):
        seq=self.t.ref[:19]+'A'+self.t.ref[19:];c='19M1D'+str(len(self.t.ref)-19)+'M'
        p=parse_record(samline(self.t,cigar=c,nm=1,score=-8),self.t,'R',self.ref(seq),CFG);self.assertEqual(p.reason,'GAP_NEAR_FOCAL_SNP')
    def test_nonallelic_reference_blocks(self):
        seq=self.t.ref[:20]+'T'+self.t.ref[21:]
        p=parse_record(samline(self.t,nm=1,score=-6),self.t,'R',self.ref(seq),CFG);self.assertEqual(p.reason,'REFERENCE_ALLELE_UNSUPPORTED')
    def test_ambiguous_reference_blocks(self):
        seq=self.t.ref[:4]+'N'+self.t.ref[5:]
        p=parse_record(samline(self.t,nm=1,score=-1),self.t,'R',self.ref(seq),CFG);self.assertEqual(p.reason,'REFERENCE_ALLELE_UNSUPPORTED')
    def test_soft_clip_blocks(self):
        c='3S'+str(len(self.t.ref)-3)+'M';p=parse_record(samline(self.t,pos=14,cigar=c),self.t,'R',self.ref(),CFG);self.assertEqual(p.reason,'NOT_FULL_END_TO_END')
    def test_bad_cigar(self):
        for c in ['0M','10Z','10Mbad','M','']:
            with self.assertRaises(DataError):cigar_ops(c)
    def test_missing_nm(self):
        s=samline(self.t).replace('\tNM:i:0','')
        with self.assertRaises(DataError):parse_record(s,self.t,'R',self.ref(),CFG)

class Decisions(unittest.TestCase):
    def test_no_alignment(self):self.assertEqual(allele_decision([],CFG)[0]['status'],'UNMAPPED')
    def test_single(self):self.assertIsNotNone(allele_decision([phit()],CFG)[1])
    def test_tie(self):self.assertEqual(allele_decision([phit(),phit(pos=100,snp=110)],CFG)[0]['status'],'INSUFFICIENT_SCORE_SEPARATION')
    def test_good_margin(self):self.assertIsNotNone(allele_decision([phit(),phit(-6,pos=100,snp=110)],CFG)[1])
    def test_cutoff_censoring(self):self.assertIsNone(allele_decision([phit(-28)],CFG)[1])
    def test_small_margin(self):self.assertIsNone(allele_decision([phit(),phit(-5,pos=100,snp=110)],CFG)[1])
    def test_cap(self):self.assertEqual(allele_decision([phit(-i*6,pos=i*100,snp=i*100+20) for i in range(10)],CFG)[0]['status'],'REPORTING_CAP_REACHED')
    def test_invalid_competitor_not_removed(self):self.assertIsNone(allele_decision([phit(),phit(pos=100,snp=None,reason='FOCAL_SNP_NOT_ON_REFERENCE_BASE')],CFG)[1])
    def test_invalid_best_not_replaced(self):self.assertIsNone(allele_decision([phit(reason='TOO_MANY_NONFOCAL_EDITS'),phit(-6,pos=100,snp=110)],CFG)[1])
    def test_xs_honored(self):self.assertIsNone(allele_decision([phit(xs=-3)],CFG)[1])
    def test_same_focal_locus_not_two_genes(self):self.assertIsNotNone(allele_decision([phit(),replace(phit(),cigar='10M1D40M',end=60)],CFG)[1])
    def test_concordant_pair(self):
        tag=Tag('q',2,'k','mk','A'*50,'A'*20+'C'+'A'*29,20)
        h=phit();d={'status':'RETAINED_BOUNDED_SEARCH'};r=combine_alleles(tag,d,h,d,h,{'I':'1'});self.assertTrue(r['accepted_for_mapping'])
    def test_discordant_pair(self):
        tag=Tag('q',2,'k','mk','A'*50,'A'*20+'C'+'A'*29,20)
        d={'status':'RETAINED_BOUNDED_SEARCH'};r=combine_alleles(tag,d,phit(),d,phit(snp=200),{'I':'1'});self.assertEqual(r['placement_status'],'ALLELE_PLACEMENTS_DISAGREE')
    def test_one_allele_insufficient(self):
        tag=Tag('q',2,'k','mk','A'*50,'A'*20+'C'+'A'*29,20)
        r=combine_alleles(tag,{'status':'UNMAPPED'},None,{'status':'RETAINED_BOUNDED_SEARCH'},phit(),{'I':'1'});self.assertFalse(r['accepted_for_mapping'])
    def test_unplaced_not_primary(self):
        tag=Tag('q',2,'k','mk','A'*50,'A'*20+'C'+'A'*29,20)
        d={'status':'RETAINED_BOUNDED_SEARCH'};r=combine_alleles(tag,d,phit(),d,phit(),{'I':None});self.assertFalse(r['accepted_for_mapping']);self.assertEqual(r['placement_status'],'RETAINED_UNPLACED_CONTIG')

class SamStream(unittest.TestCase):
    def setUp(self):
        self.t=Tag('q',2,'k','mk','ACGTTGCACTTAGCGACATCGATGTACGATTCGGAATCTAGCTACGTAG','ACGTTGCACTTAGCGACATCAATGTACGATTCGGAATCTAGCTACGTAG',20)
        self.ref=MemoryRef({'I':'N'*10+self.t.ref+'N'*20})
    def run_sam(self,s):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.sam';p.write_text(s, encoding='utf-8');return parse_sam(p,{'q':self.t},self.ref,CFG,{'I':'1'})
    def text(self):return '@SQ\tSN:I\tLN:'+str(self.ref.entries['I']['length'])+'\n'+samline(self.t)+samline(self.t,'A',nm=1,score=-6)
    def test_roundtrip(self):self.assertTrue(self.run_sam(self.text())[0][0]['accepted_for_mapping'])
    def test_missing_query_error(self):
        with self.assertRaises(DataError):self.run_sam(self.text().rsplit('\n',2)[0]+'\n')
    def test_noncontiguous_error(self):
        with self.assertRaises(DataError):self.run_sam(self.text()+samline(self.t))
    def test_wrong_dictionary(self):
        with self.assertRaises(DataError):self.run_sam(self.text().replace('LN:'+str(self.ref.entries['I']['length']),'LN:777'))
    def test_unknown_id(self):
        with self.assertRaises(DataError):self.run_sam(self.text().replace('q__R','x__R'))
