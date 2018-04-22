import numpy as np
import time
import matplotlib.pyplot as plt

def Theory(w,L):
##optimized spectral radius 
  u=np.cos(np.pi/(L-1.))
  wopt=2./(1.+np.sqrt(1.-u**2))
  if wopt <= w:
    pre_rate=w-1.
  else:
    pre_rate=1.-w+w**2*u**2/2.+w*u*np.sqrt(1-w+(w**2*u**2)/4.)
  return pre_rate

def Exact(x,y,L):
  p=np.zeros(L)
  p[:]=1
  phi=np.outer(x,p)+np.outer(x,y)+np.outer(p,y)+1
  return phi

def f(x,y):
  fxy=np.zeros((L,L))
  return fxy

def ini(L):
  x=np.linspace(0,1,L)
  y=np.linspace(0,1,L)
  ui=np.zeros((L,L))
  ui[:,0]=x+1
  ui[:,-1]=2*x+2
  ui[0,:]=y+1
  ui[-1,:]=2*y+2
  ui[1:-1,1:-1]=1
  return ui

def Jacobi(u,fxy):
	u_old[:,:]=u[:,:]
	u[1:-1, 1:-1] = (u[0:-2, 1:-1] + u[2:, 1:-1] + u[1:-1,0:-2] + u[1:-1, 2:]+h**2*fxy[1:-1,1:-1])/4.
	err=np.sum(np.sqrt((u-u_old)**2))
	return u,err

def GS_SOR_test(ux,L,w,step,stab):
  h=1./(L-1)
  x=np.linspace(0,1,L)
  y=np.linspace(0,1,L)
  fxy=f(x,y)
  err=np.zeros(step-stab)
  nstep=np.linspace(stab+1,step,num=step-stab)
  u=np.copy(ux)
  for ii in range(0,step+1):
    m, n = ux.shape
    h2 = h * h;
    for sweep in ('red', 'black'):
      if sweep == 'red':
        u[1:-1:2,1:-1:2] = (u[0:-2:2,1:-1:2] + u[2: :2,1:-1:2] + u[1:-1:2,0:-2:2] + u[1:-1:2,2: :2])/4.
        u[2:-2:2,2:-2:2] = (u[1:-3:2,2:-2:2] + u[3:-1:2,2:-2:2] + u[2:-2:2,1:-3:2] + u[2:-2:2, 3:-1:2])/4.
      else:
        u[2:-2:2,1:-1:2] = (u[1:-3:2,1:-1:2] + u[3:-1:2,1:-1:2] + u[2:-2:2,0:-2:2] + u[2:-2:2,1:-1:2])/4.
        u[1:-1:2,2:-2:2] = (u[0:-2:2,2:-2:2] + u[2: :2,2:-2:2] + u[1:-1:2,1:-3:2] + u[1:-1:2,2:-2:2])/4.  
    err[ii-stab-1]=np.sum(np.abs(u-Exact(x,y,L)))/(L-2)**2 
  fit=np.polyfit(nstep,np.log(err),2)
  rate=np.exp(fit[1])
  del ux,u
  return rate, w

def pure_python_red_black_gauss_seidel_step(u, f, h):
    m, n = u.shape
    h2 = h * h;
    for sweep in ('red', 'black'):
        for i in range(1, m - 1):
            start = 1 + i % 2 if sweep == 'red' else 2 - i % 2
            for j in range(start, n - 1, 2):
                u[i, j] = (u[i + 1, j] +
                           u[i - 1, j] +
                           u[i, j + 1] +
                           u[i, j - 1] +
                           h2 * f[i, j]) * 0.25
    return u

def GS_SOR(ux, L, w, step,stab):
  h=1./(L-1)
  x=np.linspace(0,1,L)
  y=np.linspace(0,1,L)
  fxy=f(x,y)
  err=np.zeros(step-stab)
  nstep=np.linspace(stab+1,step,num=step-stab)
  for ii in range(0,step+1):
    m, n = ux.shape
    h2 = h * h;
    for sweep in ('red', 'black'):
      for i in range(1, m - 1):
        start = 1 + i % 2 if sweep == 'black' else 2 - i % 2
        for j in range(start, n - 1, 2):
            ux[i, j] = (1-w)*ux[i,j] + w*(ux[i + 1, j] +
                      ux[i - 1, j] +
                      ux[i, j + 1] +
                      ux[i, j - 1] +
                      h2 * fxy[i, j]) * 0.25
    if ii>stab:
      err[ii-stab-1]=np.sum(np.abs(ux-Exact(x,y,L)))/(L-2)**2 

  fit=np.polyfit(nstep,np.log(err),2)
  rate=np.exp(fit[1])
  del ux
  return rate, w

#main 
start = time.time()

#Grid
step=100
stab=70
L=66
xmax=1.
ymax=1.
h=1./(L-1)
x=np.linspace(0,1,L)
y=np.linspace(0,1,L)
#Poisson function
fxy=np.zeros((L,L))
rate=np.zeros(101)
ww=np.zeros(101)
thy=np.zeros(101)

start = time.time()
#Iteration
for jj in range(0,101):
  ui=ini(L)
  w=jj/100.+1.
  result=GS_SOR(ui,L,w,step,stab)
  rate[jj]=result[0]
  ww[jj]=w
  thy[jj]=Theory(w,L)
  if jj%10==1:
    print jj


for jj in range(0,101):
  if rate[jj]==np.min(rate):
    w_opt=ww[jj]
  if thy[jj]==np.min(thy):
    w_opt_thy=ww[jj]


explot, = plt.plot(ww,rate,"r+")
thyplot, = plt.plot(ww,thy,"b-")
plt.legend([explot, thyplot],['N={:2d}'.format((L-2)),'Theory'],loc=4)
plt.axis([1,2,0.8,1.0])
plt.title("Convergence Rate v.s. $w$\n $\phi=x+x*y+y+1$")
plt.xlabel("$w$")
plt.ylabel("Convergence Rate")
elapsed = (time.time() - start)
plt.text(1.1,0.85,"time={:.3f}\nstable step={:d}\ntotal step={:d}\nw_opt={:.2f}\nw_opt_theory={:.2f}".format(elapsed,stab,step,w_opt,w_opt_thy))
plt.savefig("Conv_rate_N_{:d}_morestep.png".format(L-2))
plt.show()









