"""Cross-export identity, missing-anchor and sample-order boundary tests."""
import unittest
from types import SimpleNamespace
import numpy as np
from cross_peak_audit import geometry,classification,prepare_shared,raw_pair_correlations,best_partner
from genotype_audit import Panel

def marker(key,position,chrom='1'):
    return dict(marker_key=key,source_tag_key=key,chromosome=chrom,position_bp=position)
def peak(position):
    return dict(query_id='p'+str(position),chromosome='1',start_bp=position,end_bp=position,n_events=1,traits='test')
def physical(name,position):
    return dict(resource=name,in_C=True,reference_chromosome='1',reference_snp_bp=position)

class CrossPeakTests(unittest.TestCase):
    def test_raw_pair_formula_missing_rejection_and_variance(self):
        x=np.array([[0,1,2],[1,1,1],[0,1,2]],float);y=np.array([[2,1,0],[0,1,2],[0,2,1]],float)
        r=raw_pair_correlations(x,y,minimum=3)
        self.assertEqual(r['status'].tolist(),[0,2,0]);self.assertAlmostEqual(r['r2'][0],1.)
        self.assertAlmostEqual(r['r2'][2],np.corrcoef(x[2],y[2])[0,1]**2)
        self.assertTrue(np.isnan(r['r2'][1]))
        x[0,0]=np.nan
        with self.assertRaises(ValueError):raw_pair_correlations(x,y,minimum=3)
    def test_near_tie_uses_position_then_key(self):
        pairs=[dict(partner_position_bp=200,partner_marker_key='z'),dict(partner_position_bp=100,partner_marker_key='b'),dict(partner_position_bp=100,partner_marker_key='a')]
        r=np.array([.9,.9-2e-15,.9-3e-15])
        self.assertEqual(best_partner([0,1,2],pairs,r),2)
        self.assertIsNone(best_partner([],pairs,r))
    def test_target_anchor_precedes_donor(self):
        panels={'Nacional_G7':SimpleNamespace(markers=[marker('1|F|0--1',100),marker('2|F|0--1',150)]),
                'Nacional_MAF5':SimpleNamespace(markers=[marker('3|F|0--1',100),marker('4|F|0--1',150)])}
        rows,pairs=geometry([peak(100)],panels,[physical(n,100) for n in panels])
        self.assertEqual(rows[0]['coordinate_anchor_source'],'Nacional_G7')
        self.assertEqual(rows[1]['coordinate_anchor_source'],'Nacional_MAF5')
        self.assertEqual(rows[0]['coordinate_anchor_marker_key'],'1|F|0--1')
        self.assertTrue(all(p['distance_bp']>0 for p in pairs))
    def test_donor_boundary_and_same_key_exclusion(self):
        panels={'Nacional_G7':SimpleNamespace(markers=[marker('7|F|0--1',99999),marker('8|F|0--1',150000),marker('9|F|0--1',150001)]),
                'Nacional_MAF5':SimpleNamespace(markers=[marker('7|F|0--1',100000)])}
        rows,pairs=geometry([peak(100000)],panels,[physical('Nacional_MAF5',100000),physical('Nacional_G7',150000)])
        row=rows[0]; pp=[p for p in pairs if p['target_resource']=='Nacional_G7']
        self.assertEqual(row['anchor_selection'],'OTHER_EXPORT_QC_ANCHOR_TARGET_PHYSICALLY_ABSENT')
        self.assertEqual([p['partner_marker_key'] for p in pp],['8|F|0--1'])
        self.assertEqual(pp[0]['distance_bp'],50000)
    def test_absence_and_variance_are_not_low_correlation(self):
        r=dict(target_resource='Nacional_G7',coordinate_anchor_source='',target_shared_QC_anchor=False,target_exact_C_records=0)
        self.assertIn('NOT_ASSESSABLE',classification(r,2,0,None))
        r['coordinate_anchor_source']='Nacional_MAF5'
        self.assertIn('NO_ELIGIBLE',classification(r,0,0,None))
        self.assertIn('NOT_ASSESSABLE',classification(r,1,0,None))
        self.assertIn('BELOW_R2',classification(r,1,1,.7999))
        self.assertIn('HIGH_NEIGHBOR',classification(r,1,1,.8))
    def test_qc_excluded_exact_coordinate_remains_distinct(self):
        r=dict(target_resource='Nacional_G7',coordinate_anchor_source='Nacional_MAF5',target_shared_QC_anchor=False,target_exact_C_records=1)
        self.assertIn('QC_EXCLUDED',classification(r,1,1,.9))
        self.assertNotIn('TARGET_ABSENT',classification(r,1,1,.9))
    def test_shared_order_and_extra_sample_maf_qc(self):
        shared=['s%03d'%i for i in range(150)]+['SNA604','L17H53']
        counts={s:i%3 for i,s in enumerate(shared)}
        sa=list(reversed(shared));sb=shared+['EXTRA']
        da=np.array([[counts[s] for s in sa]],dtype=float)
        rare={s:(2 if s in shared[:7] else 0) for s in shared};rare['EXTRA']=2
        db=np.array([[counts.get(s,1) for s in sb],[rare[s] for s in sb]],dtype=float)
        pa=Panel('Nacional_G7',sa,[marker('1|F|0--1',100)],da,{'all_qc_samples':np.arange(152)}, {})
        pb=Panel('Nacional_MAF5',sb,[marker('1|F|0--1',100),marker('2|F|0--1',200)],db,{'all_qc_samples':np.arange(153)}, {})
        class Context:
            baseline_config={'dosage_ld':{'resources':['Nacional_G7','Nacional_MAF5'],'sensitivity_sample_ids':['SNA604','L17H53'],
                'qc':{'marker_min_call_rate':.9,'min_maf':.05,'min_pair_samples':50,'min_markers':1}}}
            genotype_panels={'Nacional_G7':pa,'Nacional_MAF5':pb}
            def write_csv(self,*args):pass
            def check(self,name,observed,expected):
                if observed!=expected:raise AssertionError(name)
        out=prepare_shared(Context())
        a,b=out['Nacional_G7'],out['Nacional_MAF5']
        self.assertEqual(a.samples,sorted(shared));self.assertEqual(a.samples,b.samples)
        self.assertEqual(len(b.markers),1)
        self.assertTrue(np.array_equal(a.dosage,b.dosage))
        self.assertEqual(len(a.groups['shared_omit_flagged']),150)

if __name__=='__main__':unittest.main()
