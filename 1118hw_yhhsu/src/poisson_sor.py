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

def phi(N):
  exact_2=np.zeros((N,N))
  h=1./(N-1)
  for i in range(0,N):
    for j in range(0,N):
      exact_2[i,j]=np.sin(2*np.pi*i*h)*np.sin(2*np.pi*j*h)  
  return exact_2

def f(N):
  return -8*np.pi**2*phi(N)

def ini(N):
  ui=np.zeros((N,N))
  ui[:,0]=0
  ui[:,-1]=0
  ui[0,:]=0
  ui[-1,:]=0
  ui[1:-1,1:-1]=0.5
  return ui

def GS_SOR(ux, fxy,N, w,step,stab):
  ux[:,0]=ux[:,-2]
  ux[:,-1]=ux[:,1]
  ux[0,:]=ux[-2,:]
  ux[-1,:]=ux[1,:]
  h=1./(N-1)
  for ii in range(0,1):
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
                      h2*fxy[i,j] ) * 0.25
  return ux

def residual(ux,fxy,L):
  N=2**L+2
  h=1./(N-1)
  h2=h*h
  r=np.zeros((N,N))
  r[1:-1, 1:-1] = -(ux[0:-2, 1:-1] + ux[2:, 1:-1] + ux[1:-1,0:-2] + ux[1:-1, 2:]-4*ux[1:-1,1:-1])/h2+fxy[1:-1,1:-1]
  return r
#main 

#Grid
start = time.time()
step=6
stab=0

L=5
N=2**L+2
h=1./(L-1)
w=1
#Poisson function
ux=ini(N)-0.5
fxy=f(N)
m=[]
err_plot=[]

#Iteration
for jj in range(0,10000):
  ux=GS_SOR(ux,fxy, N,w,step,stab)
  r=residual(ux,fxy,L)
  err=1./(N-2)**2*np.sum(np.abs(r))
  m.append(jj)
  err_plot.append(err)
  print err
  if err < 10**-8:
    print "err is {:.13f}".format(err)
    break
  if jj>2000:
    break

##Asymptatic fitting
elapsed = (time.time() - start)
err=np.array(err_plot)
p=np.array(m)
err=np.log(err[0::20])
n=p[0::20]
fit=np.polyfit(n[5:],err[5:],1)
slope=fit[0]
const=fit[1]
xx=np.linspace(0,1000,1000)
fitcurve=slope*xx+const

# plt.plot(err,"r+")
# plt.show()

#Plot
converge, = plt.plot(n,err,"r+")
fitting, = plt.plot(fitcurve,"b-")
plt.legend([converge, fitting],['Experiment','Fitting'],loc=1)
plt.axis([0,600,-30,5])
plt.title("Convergence rate of pure Gauss Seidel method")
plt.xlabel("N circles")
plt.ylabel("Log(Error)")
plt.text(50,-10,"time={:.3f}\nLinear Fitting\ny={:.5f} x + {:.5f}".format(elapsed,slope,const))
# plt.show()
plt.savefig("GS_converge_rate.png")

# phi_m = np.linspace(0, 1, N)
# phi_p = np.linspace(0, 1, N)
# X,Y = np.meshgrid(phi_p, phi_m)
# Z = periodic_func(X, Y).T
# exact=phi(N)

# fig = plt.figure(figsize=(14,6))

# # `ax` is a 3D-aware axis instance because of the projection='3d' keyword argument to add_subplot
# ax = fig.add_subplot(1, 2, 1, projection='3d')

# p = ax.plot_surface(X, Y, exact, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# cb = fig.colorbar(p, shrink=0.5)

# # surface_plot with color grading and color bar
# ax = fig.add_subplot(1, 2, 2, projection='3d')
# p = ax.plot_surface(X, Y, ux, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# cb = fig.colorbar(p, shrink=0.5)

# plt.show()
