"""Boundary and independent-reference checks for new CODE2 LD diagnostics."""
import unittest
from types import SimpleNamespace
import numpy as np
from ld_audit import correlation,local_pairs,neighbors,family

class LDTests(unittest.TestCase):
    def test_pair_census_against_brute_force(self):
        markers=[dict(chromosome=c,position_bp=p) for c,p in [('2',300),('1',1101),('1',100),('1',101),('2',301),('1',1100)]]
        expected=[(i,j) for i in range(len(markers)) for j in range(i+1,len(markers)) if markers[i]['chromosome']==markers[j]['chromosome'] and abs(markers[i]['position_bp']-markers[j]['position_bp'])<=1000]
        self.assertEqual(local_pairs(markers,1000).tolist(),[list(x) for x in expected])
        with self.assertRaises(ValueError):local_pairs(markers,1000,1)

    def test_scalar_reference_missing_and_allele_flip(self):
        rng=np.random.default_rng(17);d=rng.integers(0,3,(8,61)).astype(float);d[rng.random(d.shape)<.2]=np.nan
        p=np.array([(i,j) for i in range(8) for j in range(i+1,8)]);result=correlation(d,p,10,3)
        for z,(i,j) in enumerate(p):
            ok=np.isfinite(d[i])&np.isfinite(d[j]);self.assertEqual(result['n'][z],sum(ok))
            self.assertAlmostEqual(result['r2'][z],np.corrcoef(d[i,ok],d[j,ok])[0,1]**2,places=13)
        d[::2]=2-d[::2];np.testing.assert_allclose(result['r2'],correlation(d,p,10,8)['r2'],atol=1e-13,rtol=0)

    def test_missing_and_constant_are_not_low_ld(self):
        d=np.array([[0,1,np.nan,2],[1,1,1,1],[0,1,2,2.]])
        v=correlation(d,np.array([[0,1],[0,2],[1,2]]),4)
        self.assertEqual(v['status'].tolist(),[1,1,2]);self.assertTrue(np.isnan(v['r2']).all())

    def test_neighbors_no_self_stable_ties_and_exclusion(self):
        m=[dict(marker_key=str(i),chromosome='1',position_bp=p) for i,p in enumerate([100,200,300])]
        panel=SimpleNamespace(markers=m);p=np.array([[0,1],[0,2],[1,2]]);d=np.array([100,200,100])
        values=dict(r2=np.array([.8,.8,.5]),n=np.array([60]*3),status=np.array([0]*3))
        rows,summary=neighbors(panel,p,d,values,np.array([True,True,True]))
        self.assertEqual(rows[0]['best_partner_key'],'1');self.assertEqual(summary['markers_best_r2_ge_0_8'],3)
        rows,_=neighbors(panel,p,d,values,np.array([False,False,True]))
        self.assertIsNone(rows[0]['best_r2']);self.assertEqual(rows[0]['status'],'NO_OTHER_MARKER_IN_WINDOW')

    def test_source_tag_family_removes_only_terminal_index(self):
        self.assertEqual(family('3224616|F|0--36'),'3224616|F|0')
        self.assertIsNone(family('unknown'))
        self.assertEqual(family('3224616|F|0--36'),family('3224616|F|0--42'))

if __name__=='__main__':unittest.main()
