"""Tests for scientific boundary cases rather than fixed implementation counts."""
import unittest
import json
from pathlib import Path
from collections import Counter
import numpy as np
from background_audit import (gene_density, quantile_cuts, nearest_distances,
                              sample_strata, rng_for, empirical_percentile,
                              assignments, model_specs, build_catalog, SETTINGS)


def gene(gid, chrom, start, end=None, candidate=False, biotype='protein_coding'):
    end = start if end is None else end
    return dict(gene_id=gid, query_id=gid, chromosome=chrom, start_bp=start, end_bp=end,
                gene_length_bp=end-start+1, excluded_candidate=candidate, biotype=biotype)


class BackgroundTests(unittest.TestCase):
    def test_recorded_settings_match_executable_grid(self):
        saved=json.loads((Path(__file__).resolve().parent / 'SENSITIVITY_SETTINGS.json').read_text())
        self.assertEqual(saved,dict(SETTINGS,models=model_specs()))

    def test_density_boundaries_duplicate_midpoints_and_chromosomes(self):
        genes = [gene('a','1',100), gene('b','1',200), gene('c','1',200),
                 gene('d','1',301), gene('e','2',200)]
        np.testing.assert_array_equal(gene_density(genes,100), [2,2,2,0,0])
        np.testing.assert_array_equal(gene_density(genes,0), [0,1,1,0,0])

    def test_nearest_closed_interval_flanks_missing_chromosome(self):
        genes = [gene('a','1',10,20), gene('b','1',21,25), gene('c','1',1,3),gene('d','2',10,20)]
        d = nearest_distances(genes, {'1':np.array([5,20,30])})
        np.testing.assert_array_equal(d[:3], [0,1,2])
        self.assertTrue(np.isinf(d[3]))
        self.assertEqual(int(np.sum(d<=1)),2)
        with self.assertRaises(ValueError):
            nearest_distances(genes, {'1':np.array([20,5])})

    def test_tied_cutpoints_and_right_assignment(self):
        cuts = quantile_cuts([1,1,1,1],5)
        np.testing.assert_array_equal(cuts,[1])
        np.testing.assert_array_equal(np.searchsorted(cuts,[0,1,2],side='right'),[0,1,1])
        np.testing.assert_array_equal(quantile_cuts([1,2,3],1),[])

    def test_strata_without_replacement_determinism_and_deficit(self):
        pools={('1',0):np.array([0,2,4,6]),('2',0):np.array([1,3,5])}
        req=Counter({('1',0):2,('2',0):3})
        a=sample_strata(pools,req,rng_for(1729,'test'))
        b=sample_strata(pools,req,rng_for(1729,'test'))
        np.testing.assert_array_equal(a,b)
        self.assertEqual(len(a),len(set(a)))
        self.assertEqual(set(a[2:]),{1,3,5})
        with self.assertRaises(ValueError):
            sample_strata(pools,{('2',0):4},rng_for(1729,'test'))

    def test_percentile_ties_are_explicit_and_descriptive(self):
        r=empirical_percentile([2,3,3,4],3)
        self.assertEqual((r['n_draws_below'],r['n_draws_equal'],r['n_draws_above']),(1,2,1))
        self.assertEqual((r['percentile_strict'],r['percentile_midrank'],r['percentile_inclusive']),(25,50,75))
        r=empirical_percentile([3,3,3],3)
        self.assertEqual((r['percentile_strict'],r['percentile_midrank'],r['percentile_inclusive']),(0,50,100))

    def test_coding_pool_excludes_all_candidate_types_and_sparse_strata(self):
        genes=[gene('t','1',10,20,True),gene('excluded','1',100,110,True),
               gene('pool','2',100,110),gene('other','1',1000,1010,True,'nontranslating_CDS')]
        spec=next(x for x in model_specs() if x['model_id']=='L3')
        densities={250000:gene_density(genes,250000)}
        req,pools,strata,*_=assignments(genes,np.array([0]),spec,densities)
        self.assertFalse(strata[0]['adequate'])
        self.assertEqual(strata[0]['available_noncandidate_genes'],0)
        self.assertEqual(set(i for pool in pools.values() for i in pool),{2})

    def test_draw_shared_across_contrasts_uses_same_gene_indices(self):
        matrix=np.array([[0,1],[1,0],[1,1],[0,0]],dtype=bool)
        ix=sample_strata({('chr',):np.arange(4)},{('chr',):3},rng_for(1729,'shared'))
        counts=matrix[ix].sum(axis=0)
        self.assertEqual(counts.tolist(),[sum(matrix[i,j] for i in ix) for j in range(2)])

    def test_catalog_validation_and_sorting(self):
        raw=[gene('b','1',100,110),gene('a','1',1,20),gene('bad','1',1010,1020)]
        genes,excluded=build_catalog(raw,{'1':1000},{'b'})
        self.assertEqual([g['gene_id'] for g in genes],['a','b'])
        self.assertTrue(genes[1]['excluded_candidate'])
        self.assertEqual(excluded[0]['gene_id'],'bad')
        with self.assertRaises(ValueError):
            build_catalog(raw+[raw[0]],{'1':1000},set())

if __name__=='__main__':
    unittest.main()
