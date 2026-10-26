"""Independent numerical checks for conditional covariance adjustment."""
import unittest
import numpy as np
from adjusted_math import relationship_components,metric_transform,pair_correlations,metric_pairvalues,project_dosages
class AdjustedMathTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(923511)
        self.d=self.rng.integers(0,3,(25,80)).astype(float)
        self.pairs=np.array([(i,j) for i in range(25) for j in range(i+1,25)])
    def test_identity_is_pearson_and_indexed_equals_paired(self):
        t,info=metric_transform(80,V=np.eye(80))
        got=metric_pairvalues(self.d,self.pairs,t,batch=19)
        paired=pair_correlations(self.d[self.pairs[:,0]],self.d[self.pairs[:,1]],t)
        reference=np.array([np.corrcoef(self.d[i],self.d[j])[0,1]**2 for i,j in self.pairs])
        np.testing.assert_allclose(got['r2'],reference,atol=1e-12,rtol=1e-10)
        np.testing.assert_allclose(got['r2'],paired['r2'],atol=1e-12)
        self.assertEqual(info['residual_dimension'],79)
    def test_weighted_intercept_covariates_against_direct_gls(self):
        a=self.rng.normal(size=(80,90));v=a@a.T/90+.2*np.eye(80)
        covariates=self.rng.normal(size=(80,3))
        for s in (None,covariates):
            t,_=metric_transform(80,V=v,covariates=s)
            inv=np.linalg.inv(v)
            c=np.ones((80,1)) if s is None else np.column_stack((np.ones(80),s))
            precision=inv-inv@c@np.linalg.inv(c.T@inv@c)@c.T@inv
            reference=[]
            for i,j in self.pairs:
                x,y=self.d[i],self.d[j]
                reference.append((x@precision@y)**2/((x@precision@x)*(y@precision@y)))
            got=metric_pairvalues(self.d,self.pairs,t)
            np.testing.assert_allclose(got['r2'],reference,atol=1e-11,rtol=1e-10)
    def test_pc_residual_against_covariance_schur_complement(self):
        s=self.rng.normal(size=(80,5));t,_=metric_transform(80,covariates=s)
        x,y=self.d[0],self.d[1]
        full=np.cov(np.column_stack((x,y,s)),rowvar=False)
        residual=full[:2,:2]-full[:2,2:]@np.linalg.inv(full[2:,2:])@full[2:,:2]
        expected=residual[0,1]**2/(residual[0,0]*residual[1,1])
        self.assertAlmostEqual(pair_correlations(x[None,:],y[None,:],t)['r2'][0],expected,places=12)
    def test_allele_flip_covariance_scale_sample_permutation(self):
        bg=self.rng.integers(0,3,(300,80)).astype(float);obj=relationship_components(bg)
        v=obj['K']+1e-6*np.eye(80);t,_=metric_transform(80,V=v,covariates=obj['pcs'][:,:2])
        original=metric_pairvalues(self.d,self.pairs,t)['r2']
        flip=self.d.copy();flip[::2]=2-flip[::2]
        np.testing.assert_allclose(metric_pairvalues(flip,self.pairs,t)['r2'],original,atol=1e-10)
        t2,_=metric_transform(80,V=v*3,covariates=obj['pcs'][:,:2])
        np.testing.assert_allclose(metric_pairvalues(self.d,self.pairs,t2)['r2'],original,atol=1e-9)
        perm=self.rng.permutation(80);t3,_=metric_transform(80,V=v[np.ix_(perm,perm)],covariates=obj['pcs'][perm,:2])
        np.testing.assert_allclose(metric_pairvalues(self.d[:,perm],self.pairs,t3)['r2'],original,atol=1e-9)
        bg[::3]=2-bg[::3];alternate=relationship_components(bg)
        np.testing.assert_allclose(alternate['K'],obj['K'],atol=1e-12)
    def test_zero_variance_insufficient_sample_not_zero_ld(self):
        t,_=metric_transform(80);zeros=pair_correlations(np.ones((1,80)),self.d[:1],t)
        self.assertEqual(zeros['status'][0],2);self.assertTrue(np.isnan(zeros['r2'][0]))
        t,_=metric_transform(20);short=pair_correlations(self.d[:1,:20],self.d[1:2,:20],t)
        self.assertEqual(short['status'][0],1);self.assertTrue(np.isnan(short['r2'][0]))
        s=np.arange(80,dtype=float)[:,None];t,_=metric_transform(80,covariates=s)
        collapsed=pair_correlations(s.T,self.d[:1],t)
        self.assertEqual(collapsed['status'][0],2);self.assertTrue(np.isnan(collapsed['r2'][0]))
    def test_intercept_covariance_augmentation_leaves_gls_metric_unchanged(self):
        bg=self.rng.integers(0,3,(300,80)).astype(float);obj=relationship_components(bg)
        v=obj['K']+1e-6*np.eye(80)
        for covariates in (None,obj['pcs'][:,:2]):
            a,_=metric_transform(80,V=v,covariates=covariates)
            b,_=metric_transform(80,V=v+np.ones((80,80))/80,covariates=covariates)
            np.testing.assert_allclose(metric_pairvalues(self.d,self.pairs,a)['r2'],
                                       metric_pairvalues(self.d,self.pairs,b)['r2'],atol=1e-9,rtol=1e-8)
    def test_missing_invalid_matrices_rejected(self):
        t,_=metric_transform(80);bad=self.d.copy();bad[0,0]=np.nan
        with self.assertRaises(ValueError):project_dosages(bad,t)
        with self.assertRaises(ValueError):metric_transform(80,V=-np.eye(80))
        with self.assertRaises(ValueError):metric_transform(80,covariates=np.ones((80,1)))
        centered=np.eye(80)-np.ones((80,80))/80
        with self.assertRaises(ValueError):metric_transform(80,V=centered)
        with self.assertRaises(ValueError):relationship_components(np.zeros((20,80)))
        with self.assertRaises(ValueError):metric_pairvalues(self.d,np.array([[1,0]]),t)
        with self.assertRaises(ValueError):metric_pairvalues(self.d,np.array([[0,1],[0,1]]),t)
if __name__=='__main__':unittest.main()
