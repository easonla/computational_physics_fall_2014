import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def Theory(w,L):
  u=np.cos(np.pi/(L-1.))
  wopt=2./(1.+np.sqrt(1.-u**2))
  if wopt <= w:
    pre_rate=w-1.
  else:
    pre_rate=1.-w+w**2*u**2/2.+w*u*np.sqrt(1-w+(w**2*u**2)/4.)
  return pre_rate

def Exact(L):
  x=np.linspace(0,1,L)
  y=np.linspace(0,1,L)
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
      err[ii-stab-1]=np.sum(np.abs(ux-Exact(L)))/(L-2)**2 

  fit=np.polyfit(nstep,np.log(err),2)
  rate=np.exp(fit[1])
  return ux

#main 
start = time.time()

#Grid
step=100
stab=10
L=66
h=1./(L-1)
#Poisson function
rate=np.zeros(101)
ww=np.zeros(101)
thy=np.zeros(101)

start = time.time()
#Iteration
for jj in range(0,1):
  ui=ini(L)
  w=jj/100.+1.
  result=GS_SOR(ui,L,w,step,stab)
  ux =result
  # rate[jj]=result[0]
  # ww[jj]=w
  # thy[jj]=Theory(w,L)
  # del ui
  if jj%10==1:
    print jj
  

for jj in range(0,101):
  if rate[jj]==np.min(rate):
    w_opt=ww[jj]
  if thy[jj]==np.min(thy):
    w_opt_thy=ww[jj]
# print ww
# print rate
elapsed = (time.time() - start)


#Plot
# explot, = plt.plot(ww,rate,"r+")
# thyplot, = plt.plot(ww,thy,"b-")
# plt.legend([explot, thyplot],['N={:2d}'.format((L-2)),'Theory'],loc=4)
# plt.axis([1,2,0.8,1.0])
# plt.title("Convergence Rate v.s. $w$\n $\phi = x^3 + y^3$")
# plt.xlabel("$w$")
# plt.ylabel("Convergence Rate")
# plt.text(1.1,0.85,"time={:.3f}\nstable step={:d}\ntotal step={:d}\nw_opt={:.2f}\nw_opt_theory={:.2f}".format(elapsed,stab,step,w_opt,w_opt_thy))
# plt.savefig("Conv_rate_N_{:d}_2.png".format(L-2))
# plt.show()


phi_m = np.linspace(0, 1, L)
phi_p = np.linspace(0, 1, L)
X,Y = np.meshgrid(phi_p, phi_m)
exact=Exact(L)

fig = plt.figure(figsize=(14,6))

# `ax` is a 3D-aware axis instance because of the projection='3d' keyword argument to add_subplot
ax = fig.add_subplot(1, 2, 1, projection='3d')

p = ax.plot_surface(X, Y, exact, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
cb = fig.colorbar(p, shrink=0.5)

# surface_plot with color grading and color bar
ax = fig.add_subplot(1, 2, 2, projection='3d')
p = ax.plot_surface(X, Y, ux, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
cb = fig.colorbar(p, shrink=0.5)

plt.show()




