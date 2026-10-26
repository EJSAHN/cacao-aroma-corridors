"""Independent checks of prespecified training-set construction."""
import unittest
from model_setup import select_background
class BackgroundTests(unittest.TestCase):
    def test_loco_bins_and_families(self):
        markers=[dict(chromosome='1',position_bp=1,marker_key='a',source_tag_key='1|F|0--1'),
                 dict(chromosome='2',position_bp=1,marker_key='b',source_tag_key='2|F|0--1'),
                 dict(chromosome='2',position_bp=99999,marker_key='c',source_tag_key='3|F|0--1'),
                 dict(chromosome='2',position_bp=100001,marker_key='d',source_tag_key='2|F|0--5'),
                 dict(chromosome='2',position_bp=100002,marker_key='e',source_tag_key='4|F|0--5'),
                 dict(chromosome='3',position_bp=100000,marker_key='f',source_tag_key='5|F|0--1')]
        chosen=select_background(markers,'1',minimum=3)
        self.assertEqual([markers[i]['marker_key'] for i in chosen],['b','e','f'])
        rev=list(reversed(markers));chosen=select_background(rev,'1',minimum=3)
        self.assertEqual([rev[i]['marker_key'] for i in chosen],['b','e','f'])
        with self.assertRaises(ValueError):select_background(markers,'1',minimum=4)
    def test_bad_tag_family_rejected(self):
        with self.assertRaises(ValueError):select_background([dict(chromosome='2',position_bp=1,marker_key='a',source_tag_key='bad')],'1',minimum=1)
if __name__=='__main__':unittest.main()
