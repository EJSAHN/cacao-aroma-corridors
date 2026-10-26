"""Boundary checks for neighbor selection and comparability denominators."""
import unittest
from types import SimpleNamespace
import numpy as np
from within_adjusted_audit import neighbor_rows,comparison_statistics
from model_setup import select_background

class WithinTests(unittest.TestCase):
    def test_near_tie_selection_and_unavailable_are_explicit(self):
        panel=SimpleNamespace(markers=[dict(marker_key=str(i),chromosome='1',position_bp=i*100) for i in range(3)])
        pairs=np.array([[0,1],[0,2]]);dist=np.array([100,200])
        vals=dict(r2=np.array([.8,.8+2e-14]),n=np.array([60,60]),status=np.array([0,0]))
        rows,summary,best=neighbor_rows(panel,pairs,dist,vals,np.ones(2,bool),tie_tolerance=1e-12)
        self.assertEqual(rows[0]['best_partner_key'],'1')
        rows,summary,best=neighbor_rows(panel,pairs,dist,vals,np.zeros(2,bool))
        self.assertTrue(all(r['best_r2'] is None and r['status']=='NO_OTHER_MARKER_IN_WINDOW' for r in rows))
    def test_transitions_do_not_turn_missing_into_low(self):
        a=np.array([.9,.9,.2,np.nan,.9]);b=np.array([.9,.1,.9,.9,np.nan])
        r=comparison_statistics(a,b,np.ones(5,bool))
        self.assertEqual(r['valid_both'],3);self.assertEqual(r['lost_high'],1)
        self.assertEqual(r['gained_high'],1);self.assertEqual(r['raw_only_evaluable'],1)
        self.assertEqual(r['adjusted_only_evaluable'],1)
    def test_background_excludes_chr_and_repeated_bins_or_tag_families(self):
        markers=[dict(marker_key='a',chromosome='1',position_bp=1,source_tag_key='a--1'),
                 dict(marker_key='b',chromosome='2',position_bp=1,source_tag_key='b--1'),
                 dict(marker_key='c',chromosome='2',position_bp=2,source_tag_key='c--1'),
                 dict(marker_key='d',chromosome='2',position_bp=100001,source_tag_key='b--2'),
                 dict(marker_key='e',chromosome='3',position_bp=1,source_tag_key='e--1')]
        self.assertEqual(select_background(markers,'1',minimum=2),[1,4])
        with self.assertRaises(ValueError):select_background(markers,'1',minimum=3)
if __name__=='__main__':unittest.main()
