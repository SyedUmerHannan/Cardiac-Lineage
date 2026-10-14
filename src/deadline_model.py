"""Sharp finite-support deadline bounds. Synthetic/exact marginals, not fitted data."""
import numpy as np
from scipy.optimize import linprog


def inputs(a,b,H):
    a,b,H=np.asarray(a,float),np.asarray(b,float),np.asarray(H,float)
    if a.ndim!=1 or b.ndim!=1 or not len(a) or not len(b): raise ValueError('Nonempty 1D marginals required')
    if any(np.any(~np.isfinite(v)) or np.any(v<0) for v in (a,b)): raise ValueError('Invalid marginal')
    if abs(a.sum()-1)>1e-12 or abs(b.sum()-1)>1e-12: raise ValueError('Marginals must sum to one')
    if H.shape!=(len(a),len(b)) or not np.all((H==0)|(H==1)): raise ValueError('H must be a matching binary matrix')
    # Only remove input summation roundoff; no silent repair of invalid masses.
    return a/a.sum(),b/b.sum(),H


def deadline_mask(division,latency,T):
    d,l=np.asarray(division,float),np.asarray(latency,float)
    if d.ndim!=1 or l.ndim!=1 or np.any(~np.isfinite(d)) or np.any(d<0) or np.any(np.isnan(l)) or np.any(l<0):
        raise ValueError('Finite nonnegative division times; nonnegative latency or +inf required')
    if not np.isfinite(T) or T<0: raise ValueError('Finite nonnegative deadline required')
    if np.any(d>T):raise ValueError('Conditional divided cohort must have division completed by T')
    return (d[:,None]+l[None,:]<=T).astype(float)


def transport_equalities(a,b):
    m,n=len(a),len(b)
    A=np.zeros((m+n,m*n))
    for i in range(m): A[i,i*n:(i+1)*n]=1
    for j in range(n): A[m+j,j::n]=1
    return A,np.r_[a,b]


def solve(c,Aeq,beq,Aub=None,bub=None,bounds=(0,None)):
    res=linprog(c,A_eq=Aeq,b_eq=beq,A_ub=Aub,b_ub=bub,bounds=bounds,method='highs',
                options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if res.status==2:return None
    if not res.success:raise RuntimeError(res.message)
    return res


def mean_bounds(a,b,H):
    a,b,H=inputs(a,b,H);A,v=transport_equalities(a,b)
    lo=solve(2*H.ravel(),A,v);hi=solve(-2*H.ravel(),A,v)
    return float(lo.fun),float(-hi.fun)


def frontier(a,b,H,mean):
    """Sharp interval for P(Y=0) at fixed divided-lineage mean; None if infeasible.
    Includes a primal witness and floating-point dual diagnostics for the lower endpoint.
    """
    a,b,H=inputs(a,b,H);m,n=H.shape
    if not np.isfinite(mean) or not 0<=mean<=2:raise ValueError('Mean must be in [0,2]')
    A,v=transport_equalities(a,b)
    Aeq=np.zeros((len(v)+1,m*n+m));Aeq[:-1,:m*n]=A;Aeq[-1,:m*n]=2*H.ravel()
    beq=np.r_[v,mean]
    Aub=np.zeros((m,m*n+m))
    for i in range(m):Aub[i,i*n:(i+1)*n]=-2*H[i];Aub[i,m*n+i]=-1
    c=np.r_[np.zeros(m*n),np.ones(m)]
    ub=np.r_[np.full(m*n,np.inf),a]
    res=solve(c,Aeq,beq,Aub,-a,bounds=list(zip(np.zeros(len(c)),ub)))
    if res is None:return None
    z=res.x[:m*n].reshape(m,n);e=(z*H).sum(axis=1)
    dual=float(res.eqlin.marginals@beq+res.ineqlin.marginals@(-a)+res.upper.marginals[m*n:]@a)
    reduced=c-Aeq.T@res.eqlin.marginals-Aub.T@res.ineqlin.marginals-res.upper.marginals
    return {'mean':float(mean),'lower':float(res.fun),'upper':float(1-mean/2),'z':z,
            'e':e,'row_zero':res.x[m*n:],'dual_value':dual,'gap':abs(float(res.fun)-dual),
            'dual_violation':float(max(0,-min(reduced))),
            'primal_residual':float(max(np.max(abs(Aeq@res.x-beq)),np.max(np.maximum(Aub@res.x+a,0))))}


def lift(z,H,zero):
    """Construct full ordered (division,L1,L2) law from row zero masses.
    Both daughter marginals equal the supplied transport; no independence assumed.
    """
    z,H,zero=np.asarray(z,float),np.asarray(H,float),np.asarray(zero,float)
    if z.ndim!=2 or H.shape!=z.shape or zero.shape!=(z.shape[0],):raise ValueError('Incompatible lift dimensions')
    if np.any(~np.isfinite(z)) or np.any(~np.isfinite(zero)) or np.any(z<0) or not np.all((H==0)|(H==1)):
        raise ValueError('Invalid lift inputs')
    a=z.sum(axis=1);e=(z*H).sum(axis=1);m,n=z.shape
    if np.any(zero<np.maximum(0,a-2*e)-32*np.finfo(float).eps) or np.any(zero>a-e+32*np.finfo(float).eps):raise ValueError('Infeasible lift')
    # Snap solver roundoff to the admissible status interval before constructing
    # the law. This changes supplied numerical marginals only at roundoff scale.
    zero=np.minimum(np.maximum(0,a-e),np.maximum(np.maximum(0,a-2*e),zero))
    full=np.zeros((m,n,n))
    for i in range(m):
        succ=z[i]*H[i];fail=z[i]*(1-H[i]);s=succ.sum();f=fail.sum()
        if s>0:succ=succ/s
        if f>0:fail=fail/f
        p00=zero[i];p11=max(0.,2*e[i]-a[i]+p00);p01=max(0.,a[i]-e[i]-p00)
        full[i]=p00*np.outer(fail,fail)+p11*np.outer(succ,succ)+p01*(np.outer(succ,fail)+np.outer(fail,succ))
        mass=full[i].sum()
        if mass>0:full[i]*=a[i]/mass
    return full


def full_oracle(a,b,H,mean,objective='min_zero'):
    """Independent uncompressed three-marginal LP used only for verification."""
    a,b,H=inputs(a,b,H);m,n=H.shape
    idx=list(np.ndindex(m,n,n));A=np.zeros((m+2*n+1,len(idx)));v=np.r_[a,b,b,mean]
    costs=[]
    for k,(i,j,l) in enumerate(idx):
        A[i,k]=1;A[m+j,k]=1;A[m+n+l,k]=1;A[-1,k]=H[i,j]+H[i,l]
        costs.append((1-H[i,j])*(1-H[i,l]))
    c=np.array(costs);sign=-1 if objective=='max_zero' else 1
    r=solve(sign*c,A,v)
    return None if r is None else float(sign*r.fun)
