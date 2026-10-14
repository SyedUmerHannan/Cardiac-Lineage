import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from deadline_model import mean_bounds,frontier,full_oracle,lift,deadline_mask

class DeadlineTests(unittest.TestCase):
 def test_closed_form_two_bin(self):
  H=np.array([[1,1],[1,0]])
  for a in np.linspace(.1,.9,9):
   for b in np.linspace(.1,.9,9):
    lo,hi=mean_bounds([a,1-a],[b,1-b],H)
    self.assertAlmostEqual(lo,2*max(a,b));self.assertAlmostEqual(hi,2*min(1,a+b))
    for mu in np.linspace(lo,hi,5):
     f=frontier([a,1-a],[b,1-b],H,mu)
     self.assertAlmostEqual(f['lower'],max(0,1+a-mu),places=9)
     self.assertAlmostEqual(f['upper'],1-mu/2,places=9)
 def test_reduced_matches_full_and_constructive_lift(self):
  rng=np.random.default_rng(20261014)
  for iteration in range(200):
   m,n=rng.integers(1,5,2);a=rng.dirichlet(np.ones(m));b=rng.dirichlet(np.ones(n))
   # Include non-threshold masks: reduction does not rely on monotone times.
   H=rng.integers(0,2,(m,n));lo,hi=mean_bounds(a,b,H)
   for mu in [lo,(lo+hi)/2,hi]:
    mu=float(np.clip(mu,0,2));f=frontier(a,b,H,mu)
    self.assertAlmostEqual(f['lower'],full_oracle(a,b,H,mu),places=8)
    self.assertAlmostEqual(f['upper'],full_oracle(a,b,H,mu,'max_zero'),places=8)
    for v in ['gap','dual_violation','primal_residual']:self.assertLess(f[v],1e-8)
    for target in [f['row_zero'],a-f['e'],.5*(f['row_zero']+a-f['e'])]:
     x=lift(f['z'],H,target)
     self.assertGreaterEqual(x.min(),0.)
     np.testing.assert_allclose(x.sum(axis=(1,2)),a,atol=1e-9)
     np.testing.assert_allclose(x.sum(axis=(0,2)),b,atol=1e-9)
     np.testing.assert_allclose(x.sum(axis=(0,1)),b,atol=1e-9)
     yy=H[:,:,None]+H[:,None,:]
     self.assertAlmostEqual(float((x*yy).sum()),mu,places=8)
     self.assertAlmostEqual(float((x*(yy==0)).sum()),float(sum(target)),places=8)
 def test_exact_boundary_and_failure_atom(self):
  H=deadline_mask([0,2],[0,2,np.inf],2)
  np.testing.assert_array_equal(H,[[1,1,0],[1,0,0]])
  self.assertEqual(mean_bounds([0,1],[0,0,1],H),(0.,0.))
 def test_faster_marginal_ranking_reversal(self):
  A=deadline_mask([0,2],[0,2],2.5);B=deadline_mask([.2,2.2],[.2,2.2],2.5)
  za=np.diag([.5,.5]);zb=np.array([[0,.5],[.5,0]])
  self.assertEqual(2*(za*A).sum(),1);self.assertEqual(2*(zb*B).sum(),2)
 def test_infeasible_and_invalid(self):
  with self.assertRaises(ValueError):deadline_mask([3],[0],2)
  with self.assertRaises(ValueError):lift([[1]],[[1]],[-1e-10])
  self.assertIsNone(frontier([.5,.5],[.5,.5],[[1,1],[1,0]],.8))
  for a,b,H in [([-.1,1.1],[1],[[1],[0]]),([.9],[1],[[1]]),([1],[1],[[.5]])]:
   with self.assertRaises(ValueError):mean_bounds(a,b,H)
 def test_paired_transport_reduces_uncertainty(self):
  z=np.array([[.25,.25],[.25,.25]]);H=np.array([[1,1],[1,0]])
  e=(z*H).sum(axis=1);a=z.sum(axis=1)
  self.assertEqual(2*e.sum(),1.5)
  self.assertEqual(np.maximum(0,a-2*e).sum(),0)
  self.assertEqual(1-e.sum(),.25)

if __name__=='__main__':unittest.main()
