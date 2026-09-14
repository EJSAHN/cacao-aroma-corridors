"""Independent numerical and pair-inventory controls, not cacao observations."""
import unittest
import numpy as np
from cacao_marker_harmonization.ld.common import DataError
from cacao_marker_harmonization.ld.numerics import pairwise_ld,near_pair_codes,sample_context_codes,pair_inventory,summary_statistics,array_digest
from cacao_marker_harmonization.ld.analysis import geometry,neighbor_diagnostics
from cacao_marker_harmonization.ld.genotypes import Panel

class CorrelationTests(unittest.TestCase):
    def test_perfect_positive(self):
        r=pairwise_ld(np.array([[0,1,2,1],[0,1,2,1.]]),np.array([[0,1]]),3)
        self.assertEqual(r.r2[0],1)
    def test_perfect_negative(self):
        r=pairwise_ld(np.array([[0,1,2,1],[2,1,0,1.]]),np.array([[0,1]]),3)
        self.assertEqual(r.r2[0],1)
    def test_zero(self):
        r=pairwise_ld(np.array([[0,0,2,2],[0,2,0,2.]]),np.array([[0,1]]),3)
        self.assertEqual(r.r2[0],0)
    def test_pairwise_missing_intersection(self):
        a=np.array([[0,1,2,np.nan,2],[0,np.nan,2,1,1.]])
        r=pairwise_ld(a,np.array([[0,1]]),3)
        self.assertEqual(r.n[0],3)
        self.assertAlmostEqual(r.r2[0],np.corrcoef(a[0,[0,2,4]],a[1,[0,2,4]])[0,1]**2)
    def test_insufficient(self):
        r=pairwise_ld(np.array([[0,1,2,np.nan],[0,1,2,2.]]),np.array([[0,1]]),4)
        self.assertEqual(r.status[0],1);self.assertTrue(np.isnan(r.r2[0]))
    def test_zero_variance(self):
        r=pairwise_ld(np.array([[1,1,1,1],[0,1,2,2.]]),np.array([[0,1]]),3)
        self.assertEqual(r.status[0],2);self.assertTrue(np.isnan(r.r2[0]))
    def test_zero_pairwise_variance(self):
        r=pairwise_ld(np.array([[0,0,0,2],[0,1,2,np.nan]]),np.array([[0,1]]),3)
        self.assertEqual(r.status[0],2)
    def test_all_missing(self):
        r=pairwise_ld(np.full((2,4),np.nan),np.array([[0,1]]),3)
        self.assertEqual(r.n[0],0);self.assertEqual(r.status[0],1)
    def test_random_against_numpy(self):
        rng=np.random.default_rng(17);d=rng.integers(0,3,(21,71)).astype(float);d[rng.random(d.shape)<.17]=np.nan
        p=np.array([(i,j) for i in range(21) for j in range(i+1,21)])
        r=pairwise_ld(d,p,20,batch=13)
        for k,(i,j) in enumerate(p):
            ok=np.isfinite(d[i])&np.isfinite(d[j])
            self.assertEqual(r.n[k],ok.sum());self.assertAlmostEqual(r.r2[k],np.corrcoef(d[i,ok],d[j,ok])[0,1]**2,places=13)
    def test_batch_independence(self):
        d=np.random.default_rng(7).integers(0,3,(10,51)).astype(float);p=np.array([(i,j) for i in range(10) for j in range(i+1,10)])
        np.testing.assert_equal(pairwise_ld(d,p,3,1).r2,pairwise_ld(d,p,3,100).r2)
    def test_allele_flip_invariant(self):
        d=np.random.default_rng(9).integers(0,3,(5,35)).astype(float);p=np.array([[0,1],[0,2],[2,4]])
        a=pairwise_ld(d,p,3).r2;d[0]=2-d[0];d[2]=2-d[2]
        np.testing.assert_allclose(a,pairwise_ld(d,p,3).r2,rtol=0,atol=1e-14)
    def test_sample_order_invariant(self):
        d=np.array([[0,1,0,2,1],[0,2,2,1,1.]])
        np.testing.assert_allclose(pairwise_ld(d, np.array([[0,1]]),3).r2,pairwise_ld(d[:,::-1],np.array([[0,1]]),3).r2)
    def test_invalid_dosage(self):
        with self.assertRaises(DataError):pairwise_ld(np.array([[0,3,1],[0,1,2.]]),np.array([[0,1]]),3)
    def test_infinite_error(self):
        with self.assertRaises(DataError):pairwise_ld(np.array([[0,np.inf,1],[0,1,2.]]),np.array([[0,1]]),3)
    def test_self_pair_rejected(self):
        with self.assertRaises(DataError):pairwise_ld(np.ones((2,4)),np.array([[1,1]]),3)
    def test_reversed_pair_rejected(self):
        with self.assertRaises(DataError):pairwise_ld(np.ones((2,4)),np.array([[1,0]]),3)
    def test_duplicate_pair_rejected(self):
        with self.assertRaises(DataError):pairwise_ld(np.ones((2,4)),np.array([[0,1],[0,1]]),3)
    def test_pair_outside_bounds(self):
        with self.assertRaises(DataError):pairwise_ld(np.ones((2,4)),np.array([[0,2]]),3)
    def test_no_pairs(self):
        self.assertEqual(len(pairwise_ld(np.ones((2,4)),np.empty((0,2),int),3).r2),0)
    def test_summary_missing_not_zero(self):
        v=pairwise_ld(np.ones((2,4)),np.array([[0,1]]),3)
        s=summary_statistics(v,np.array([True]),[.8])
        self.assertEqual(s['valid_pairs'],0);self.assertIsNone(s['mean_r2']);self.assertIsNone(s['fraction_r2_ge_0_8'])

class PairInventoryTests(unittest.TestCase):
    def test_boundary_inclusive(self):
        self.assertEqual(near_pair_codes(['1']*3,[1,11,12],10,100).tolist(),[1,5])
    def test_zero_distance_separate(self):
        self.assertEqual(near_pair_codes(['1']*2,[9,9],0,100).tolist(),[1])
    def test_different_chromosomes(self):
        self.assertEqual(len(near_pair_codes(['1','2'],[1,1],100,100)),0)
    def test_unsorted(self):
        c=['1','2','1','1'];p=[50,10,30,60];got=near_pair_codes(c,p,20,100)
        expected=[i*4+j for i in range(4) for j in range(i+1,4) if c[i]==c[j] and abs(p[i]-p[j])<=20]
        self.assertEqual(got.tolist(),expected)
    def test_random_exhaustive(self):
        rng=np.random.default_rng(11);c=rng.integers(1,4,40).astype(str);p=rng.integers(1,500,40)
        got=near_pair_codes(c,p,40,10000)
        expected=[i*40+j for i in range(40) for j in range(i+1,40) if c[i]==c[j] and abs(p[i]-p[j])<=40]
        self.assertEqual(got.tolist(),expected)
    def test_safety_limit_not_truncation(self):
        with self.assertRaises(DataError):near_pair_codes(['1']*10,list(range(1,11)),100,20)
    def test_inter_context_census(self):
        c=['1','1','2','3'];p=[1,2,3,4];got,n=sample_context_codes(c,p,100,20,1,'x','interchromosomal')
        e=[i*4+j for i in range(4) for j in range(i+1,4) if c[i]!=c[j]]
        self.assertEqual(got.tolist(),e);self.assertEqual(n,len(e))
    def test_far_context_census(self):
        c=['1','1','2','1','2'];p=[50,1,40,12,1];got,n=sample_context_codes(c,p,100,30,1,'x','distant_same_chromosome')
        e=[i*5+j for i in range(5) for j in range(i+1,5) if c[i]==c[j] and abs(p[i]-p[j])>=30]
        self.assertEqual(got.tolist(),e);self.assertEqual(n,len(e))
    def test_context_deterministic(self):
        c=['1']*30+['2']*30;p=list(range(1,61));a,_=sample_context_codes(c,p,15,10,17,'x','interchromosomal');b,_=sample_context_codes(c,p,15,10,17,'x','interchromosomal')
        np.testing.assert_equal(a,b);self.assertEqual(len(set(a)),15)
    def test_zero_context(self):
        a,n=sample_context_codes(['1','2'],[1,2],0,10,1,'x','interchromosomal');self.assertEqual(len(a),0);self.assertEqual(n,1)
    def test_empty_context(self):
        a,n=sample_context_codes(['1','1'],[1,2],10,10,1,'x','interchromosomal');self.assertEqual(len(a),0);self.assertEqual(n,0)
    def test_union_old_and_new(self):
        cfg=dict(max_distance_bp=10,max_union_pairs=100,context_pairs_per_class=0,distant_same_chromosome_min_bp=100)
        p,f,stats=pair_inventory(['1']*3,[1,2,30],['1']*3,[1,30,35],cfg,1,'x')
        self.assertEqual(p.tolist(),[[0,1],[1,2]]);self.assertEqual(stats['original_local_pairs'],1);self.assertEqual(stats['corrected_local_pairs'],1)
    def test_original_chromosome_change(self):
        same,dist=geometry(np.array([[0,1]]),np.array(['1','2']),np.array([5,10]));self.assertFalse(same[0])
    def test_digest_shape(self):
        self.assertNotEqual(array_digest(np.array([1,2])),array_digest(np.array([[1,2]])))

class NeighborTests(unittest.TestCase):
    def test_no_self_and_absent_neighbor(self):
        markers=[dict(marker_key=str(i),chromosome='1',position_bp=p) for i,p in enumerate([1,10,1000])]
        panel=Panel('x',list('abcd'),markers,np.array([[0,1,2,0],[0,1,2,0],[2,1,1,0.]]),{}, {})
        pairs=np.array([[0,1]])
        val=pairwise_ld(panel.dosage,pairs,3)
        res=dict(panel=panel,pairs=pairs,geometries={'corrected_snp':(np.array([True]),np.array([9]))},values={'all_qc_samples':val})
        rows,ss=neighbor_diagnostics(res,{'ld':{'report_windows_bp':[10],'r2_thresholds':[.8]}})
        self.assertEqual(rows[0]['best_partner_key'],'1');self.assertEqual(rows[1]['best_partner_key'],'0')
        self.assertIsNone(rows[2]['best_r2']);self.assertEqual(ss[0]['markers_with_valid_pair'],2)
        self.assertAlmostEqual(ss[0]['fraction_all_markers_ge_0_8'],2/3);self.assertEqual(ss[0]['fraction_evaluable_ge_0_8'],1)
