import unittest
import numpy as np
from cacao_marker_harmonization.coverage.common import DataError
from cacao_marker_harmonization.coverage.numerics import (Queries, positions, nearest, counts, paired_counts,
 rng_for, quotas_by_chromosome, sample_points, gene_neighborhood_counts, bin_cuts,bin_of,sample_strata)


def query(intervals):
    return Queries('test',[dict(query_id=str(i),chromosome=c,start_bp=s,end_bp=e) for i,(c,s,e) in enumerate(intervals)])


class DistanceTests(unittest.TestCase):
    def test_closed_gene_body(self):
        q=query([('1',10,20)])
        for p in [10,15,20]:self.assertEqual(nearest(q,{'1':np.array([p])})[0],0)
    def test_point_exact(self):self.assertEqual(nearest(query([('1',10,10)]),{'1':np.array([10])})[0],0)
    def test_flank_boundaries(self):
        d=nearest(query([('1',10,20)]),{'1':np.array([25])})
        np.testing.assert_array_equal(counts(d,[0,4,5,10]),[0,0,1,1])
    def test_marker_left(self):self.assertEqual(nearest(query([('1',10,20)]),{'1':np.array([7])})[0],3)
    def test_marker_right(self):self.assertEqual(nearest(query([('1',10,20)]),{'1':np.array([25])})[0],5)
    def test_missing_chromosome(self):self.assertTrue(np.isinf(nearest(query([('2',10,20)]),{'1':np.array([15])})[0]))
    def test_no_markers(self):np.testing.assert_array_equal(counts(nearest(query([('1',1,10)]),{}),[0,100]),[0,0])
    def test_empty_queries(self):self.assertEqual(len(nearest(query([]),{})),0)
    def test_duplicate_query_id(self):
        with self.assertRaises(DataError):Queries('x',[dict(query_id='a',chromosome='1',start_bp=1,end_bp=2)]*2)
    def test_invalid_interval(self):
        with self.assertRaises(DataError):query([('1',10,5)])
    def test_zero_start(self):
        with self.assertRaises(DataError):query([('1',0,5)])
    def test_unique_positions(self):
        p=positions([dict(chromosome='1',position_bp=v) for v in [3,2,3,4]])
        np.testing.assert_array_equal(p['1'],[2,3,4])
    def test_random_bruteforce(self):
        rng=np.random.default_rng(12)
        for _ in range(40):
            pp={c:np.unique(rng.integers(1,1000,size=20)) for c in ['1','2']}
            qq=[(str(rng.integers(1,4)),int(s),int(s)+int(rng.integers(0,100))) for s in rng.integers(1,1000,size=30)]
            d=nearest(query(qq),pp)
            brute=[min([max(s-int(p),int(p)-e,0) for p in pp.get(c,[])],default=np.inf) for c,s,e in qq]
            np.testing.assert_array_equal(d,brute)
    def test_monotonicity(self):
        d=np.array([0,1,5,100,np.inf]);cc=counts(d,[0,1,10,100,200]);self.assertTrue(np.all(np.diff(cc)>=0))
    def test_paired_gain_loss(self):
        s=paired_counts(np.array([0,0,0,0]),np.array([0,20,30,0]),np.array([40,0,0,0]),10)
        self.assertEqual(s['gained_after_coordinate_change'],2);self.assertEqual(s['lost_after_coordinate_change'],1)
        self.assertEqual(s['coordinate_delta_count'],1)
    def test_selection_loss_nonpositive(self):
        s=paired_counts(np.array([0,0]),np.array([0,30]),np.array([0,0]),5)
        self.assertEqual(s['selection_delta_count'],-1);self.assertEqual(s['total_delta_count'],0)
    def test_B_not_subset_error(self):
        with self.assertRaises(DataError):paired_counts(np.array([20]),np.array([0]),np.array([0]),5)
    def test_pair_shape_error(self):
        with self.assertRaises(DataError):paired_counts(np.array([0,0]),np.array([0]),np.array([0]),5)
    def test_infinite_missing_remains_miss(self):
        s=paired_counts(np.array([np.inf]),np.array([np.inf]),np.array([np.inf]),10000)
        self.assertEqual(s['neither_B_C'],1)


class SamplingTests(unittest.TestCase):
    def test_seed_stable(self):np.testing.assert_array_equal(rng_for(10,'a').integers(100,size=20),rng_for(10,'a').integers(100,size=20))
    def test_seed_streams_differ(self):self.assertFalse(np.array_equal(rng_for(10,'a').integers(100,size=20),rng_for(10,'b').integers(100,size=20)))
    def test_quotas(self):self.assertEqual(quotas_by_chromosome({'a':{'1':[1,2],'2':[3]},'b':{'1':[1]}},['1','2']),{'1':1,'2':0})
    def test_sample_without_replacement(self):
        s=sample_points({'1':np.arange(1,20)},{'1':7},rng_for(1,'x'));self.assertEqual(len(s['1']),7);self.assertEqual(len(set(s['1'])),7)
    def test_full_quota_unchanged(self):
        s=sample_points({'1':np.array([1,4,6])},{'1':3},rng_for(1,'x'));np.testing.assert_array_equal(s['1'],[1,4,6])
    def test_zero_quota(self):self.assertEqual(sample_points({}, {'1':0},rng_for(1,'a')), {})
    def test_invalid_quota(self):
        with self.assertRaises(DataError):sample_points({'1':np.array([1])},{'1':2},rng_for(1,'a'))
    def test_neighborhood_count(self):
        rr=[dict(chromosome=c,start_bp=s,end_bp=s) for c,s in [('1',10),('1',20),('1',40),('2',10)]]
        np.testing.assert_array_equal(gene_neighborhood_counts(rr,10),[1,1,0,0])
    def test_equal_midpoint_other_genes(self):
        rr=[dict(chromosome='1',start_bp=1,end_bp=5)]*3
        np.testing.assert_array_equal(gene_neighborhood_counts(rr,0),[2,2,2])
    def test_strata_exact_counts(self):
        p={('1','p',1):np.arange(10),('2','p',1):np.arange(10,20)}
        x=sample_strata(p,{('1','p',1):3,('2','p',1):4},rng_for(5,'x'))
        self.assertEqual(int((x<10).sum()),3);self.assertEqual(len(x),7)
    def test_insufficient_strata_blocks(self):
        with self.assertRaises(DataError):sample_strata({('a',):np.arange(2)},{('a',):3},rng_for(1,'x'))
    def test_quantile_bin_ties(self):
        c=bin_cuts([10,10,10,20,20],5);self.assertEqual(bin_of(10,c),bin_of(10,c));self.assertTrue(np.all(np.diff(c)>0))
    def test_empty_bins_error(self):
        with self.assertRaises(DataError):bin_cuts([],5)
    def test_one_bin(self):self.assertEqual(len(bin_cuts([1,2],1)),0)
