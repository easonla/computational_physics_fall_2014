import numpy as np
import time
import matplotlib.pyplot as plt
startup = time.time()
def Theory(w,L):
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
  phi=np.outer(x**3,p)+np.outer(p,y**3)
  return phi

def f1(x,y,L):
  p=np.zeros(L)
  p[:]=1
  fxy=np.zeros((L,L))
  fxy=6*(np.outer(x,p)+np.outer(p,y))
  return fxy


def ini(L):
  x=np.linspace(0,1,L)
  y=np.linspace(0,1,L)
  ui=np.zeros((L,L))
  ui[:,0]=x**3
  ui[:,-1]=1+x**3
  ui[0,:]=y**3
  ui[-1,:]=1+y**3
  ui[1:-1,1:-1]=1
  return ui

def GS_SOR(ux, L, w, step,stab):
  h=1./(L-1)
  x=np.linspace(0,1,L)
  y=np.linspace(0,1,L)
  fxy=f1(x,y,L)
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
                      ux[i, j - 1] -
                      6* h**3 * (i+j) ) * 0.25
    if ii>stab:
      err[ii-stab-1]=np.sum(np.abs(ux-Exact(x,y,L)))/(L-2)**2 

  fit=np.polyfit(nstep,np.log(err),2)
  rate=np.exp(fit[1])
  del ux
  return rate, w

#main 
start = time.time()

#Grid
step=50
stab=10
w_opt=np.zeros(10)
w_opt_thy=np.zeros(10)
start=np.zeros(10)
elapsed=np.zeros(10)
Grid=np.array([4,8,12,16,20,24,28,32,48,64])
for pp in range(10):
  L=Grid[pp]+2
  start[pp] = time.time()
  h=1./(L-1)
  #Poisson function
  rate=np.zeros(101)
  ww=np.zeros(101)
  thy=np.zeros(101)
  #Iteration
  for jj in range(0,101):
    ui=ini(L)
    w=jj/100.+1.
    result=GS_SOR(ui,L,w,step,stab)
    rate[jj]=result[0]
    ww[jj]=w
    thy[jj]=Theory(w,L)
    del ui

  for jj in range(0,101):
    if rate[jj]==np.min(rate):
      w_opt[pp]=ww[jj]
  elapsed[pp] = (time.time() - start[pp])

x=np.linspace(5,80,200)
u=np.cos(np.pi/(x+1.))
w_opt_thy=2./(1.+np.sqrt(1.-u**2))
endup=(time.time()-startup)
print endup

#Plot
explot, = plt.plot(Grid,w_opt,"r+")
thyplot, = plt.plot(x,w_opt_thy,"b-")
plt.legend([explot, thyplot],['optimal w','Theory'],loc=4)
plt.axis([0,70,1.49,2.0])
plt.title("optimal $w$ v.s Lattice size")
plt.xlabel("N")
plt.ylabel("$w_opt$")
plt.savefig("optimal.png".format(L-2))
plt.show()








