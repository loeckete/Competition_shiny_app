"""Independent SciPy reference and Python check of the R integrator's numerical scheme."""
import numpy as np,csv,json
from scipy.integrate import solve_ivp
results=[]
with open('tests/reference.csv','w') as f:
 w=csv.writer(f);w.writerow(['case','time','N1','N2'])
 for name,a,b,n1,n2 in [('A',.5,.5,20,20),('B',.5,1.5,20,20),('C',1.5,.5,20,20),('D',1.5,1.5,30,20),('E',1.5,1.5,20,30)]:
  def rhs(t,y):
   x,z=y;return np.array([.2*x*(1-(x+a*z)/100),.2*z*(1-(z+b*x)/100)])
  t=np.arange(0,150.01,.25)
  sol=solve_ivp(rhs,(0,150),[n1,n2],t_eval=t,method='DOP853',rtol=1e-12,atol=1e-12)
  y=np.array([n1,n2],dtype=float);out=[y.copy()]
  for _ in t[1:]:
   h=.25/3
   for j in range(3):
    k1=rhs(0,y);k2=rhs(0,y+h*k1/2);k3=rhs(0,y+h*k2/2);k4=rhs(0,y+h*k3)
    y+=h*(k1+2*k2+2*k3+k4)/6
   out.append(y.copy())
  err=float(np.max(np.abs(np.array(out)-sol.y.T)));assert err<1e-5
  assert np.min(out)>=0
  for j in [0,40,200,400,600]:w.writerow([name,t[j],*sol.y[:,j]])
  results.append(dict(case=name,maximum_error=err,final=sol.y[:,-1].tolist()))
json.dump({'scope':'Python reproduction, not native R execution','results':results},open('tests/numerical-check.json','w'),indent=2)
print(results)
