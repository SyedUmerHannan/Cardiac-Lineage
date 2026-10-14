from pathlib import Path
import csv,json,subprocess,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from deadline_model import frontier,mean_bounds,lift,deadline_mask
ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
H=np.array([[1,1],[1,0]]);a=b=np.array([.5,.5]);rows=[]
for mu in np.linspace(1,2,101):
 f=frontier(a,b,H,float(mu));rows.append({k:f[k] for k in ['mean','lower','upper','dual_value','gap','dual_violation','primal_residual']})
with (ROOT/'results/frontier.csv').open('w',newline='') as out:
 w=csv.DictWriter(out,rows[0].keys());w.writeheader();w.writerows(rows)
f=frontier(a,b,H,1.2)
for name,z0 in [('lower',f['row_zero']),('upper',a-f['e'])]:
 x=lift(f['z'],H,z0)
 with (ROOT/f'results/witness_{name}.csv').open('w',newline='') as out:
  w=csv.writer(out);w.writerow(['division_bin','daughter1_latency_bin','daughter2_latency_bin','probability','deadline_output'])
  for i,j,k in np.ndindex(x.shape):w.writerow([i,j,k,float(x[i,j,k]),int(H[i,j]+H[i,k])])
fig,ax=plt.subplots(figsize=(6.7,4.1));mu=np.array([r['mean'] for r in rows]);low=np.array([r['lower'] for r in rows]);high=np.array([r['upper'] for r in rows])
ax.fill_between(mu,low,high,color='#149c95',alpha=.28,label='Attainable with shared division time')
ax.plot(mu,low,color='#087a75',lw=2,label='Sharp lower boundary')
ax.plot(mu,high,color='#273d64',lw=2,label='Sharp upper boundary')
ax.plot(mu,np.maximum(0,1-mu),'--',color='#bd7443',label='Mean-only lower bound')
ax.plot([1.2,1.2],[.3,.4],color='black',lw=3);ax.scatter([1.2],[.1],marker='x',color='#b82b3d',s=65,label='Impossible pair (1.2, 0.1)')
ax.set(xlabel='Mean qualified daughters by deadline',ylabel='Probability neither daughter qualifies',xlim=(1,2),ylim=(-.025,.55),title='Same timing marginals constrain the joint endpoint')
ax.legend(fontsize=8);fig.savefig(ROOT/'figures/frontier.pdf');fig.savefig(ROOT/'figures/frontier.png',dpi=180);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(7,3.3))
for ax,z,title in zip(axs,[np.diag([.5,.5]),np.array([[0,.5],[.5,0]])],['A: faster marginal times; mean 1','B: slower marginal times; mean 2']):
 ax.imshow(z,cmap='Blues',vmin=0,vmax=.5)
 for i,j in np.ndindex(z.shape):ax.text(j,i,str(z[i,j]),ha='center',va='center',color='white' if z[i,j] else 'black')
 ax.set(xticks=[0,1],yticks=[0,1],xticklabels=['Fast','Slow'],yticklabels=['Early','Late'],xlabel='Rematuration latency',ylabel='Division time',title=title)
fig.tight_layout();fig.savefig(ROOT/'figures/ranking.pdf');fig.savefig(ROOT/'figures/ranking.png',dpi=180);plt.close(fig)
report={'data_status':'synthetic; no empirical calibration','two_bin_example':{'mean':1.2,'zero_probability_interval':[f['lower'],f['upper']],'mean_only_interval':[0,.4]},'frontier_grid_points':101,'max_dual_gap':max(r['gap'] for r in rows),'max_primal_residual':max(r['primal_residual'] for r in rows),'max_dual_violation':max(r['dual_violation'] for r in rows),'verification_design':{'random_masks':200,'mean_points_per_mask':3,'full_LP_comparisons':1200,'constructive_lifts':1800,'two_bin_parameter_pairs':81,'two_bin_frontier_checks':405}}
(ROOT/'results/summary.json').write_text(json.dumps(report,indent=2)+'\n')
r=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT/'tests'),'-v'],capture_output=True,text=True)
(ROOT/'results/test_log.txt').write_text(r.stdout+r.stderr)
print(json.dumps(report,indent=2));print(r.stderr)
if r.returncode:raise SystemExit(r.returncode)
