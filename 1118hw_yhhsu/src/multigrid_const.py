# Computational physics homework 6
# Multi-grid algorithm with periodic boundary condition 
# by yihsuan hsu 11/15/2014
import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
def phi(L):
##defines exact solution phi = sin (2 pi x) sin (2 pi y)
##L defines the order of grid number in each side
##there are totally 2^L unknows and an extra coloum in each side
  N=2**L+2  
  exact_2=np.zeros((N,N))
  h=1./(N-1)
  for i in range(0,N):
    for j in range(0,N):
      exact_2[i,j]=np.sin(2*np.pi*i*h)*np.sin(2*np.pi*j*h)  
  return exact_2

def f(L):
##define right hand side of poisson equation
  return -8*np.pi**2*phi(L)

def ini(L): 
##initial guess matrix is zero matrix
  N=2**L+2
  ui=np.zeros((N,N))
  # ui[1:-1,1:-1]=np.random.random((N-2, N-2))
  return ui

def GS_SOR(ux, fxy,L, w,step):
##Gauss Seidel algorithm  with Overrelaxation and periodic boundary condition
##Red-Black ordering 
##periodic BC, exchange the ghost grid to the coloum in opposite side 
  ux[:,0]=ux[:,-2]
  ux[:,-1]=ux[:,1]
  ux[0,:]=ux[-2,:]
  ux[-1,:]=ux[1,:]
  N=2**L+2
  h=1./(N-1)
  for ii in range(0,step):
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

def restrict(A,L):
##restrict L->L-1
  n=2**(L-1)+2
  AA=np.zeros((n,n))
  for i in range(1,n-1):
    for j in range(1,n-1):
      c,d=2*i-1,2*j-1
      AA[i,j]=0.25*(A[c,d]+A[c+1,d]+A[c,d+1]+A[c+1,d+1])
  return AA

def relax(A,L):
##relax L->L+1
  n=2**(L+1)+2
  AA=np.zeros((n,n))
  for i in range(1,n-1,2):
    for j in range(1,n-1,2):
      AA[i,j]=A[i/2+1,j/2+1]
      AA[i+1,j]=A[i/2+1,j/2+1]
      AA[i,j+1]=A[i/2+1,j/2+1]
      AA[i+1,j+1]=A[i/2+1,j/2+1]
  return AA


def mgm(u,f,L,w,pre_step,post_step):
##Multi-grid method

  # print"now is running L={:d}".format(L)
  # print"dimension check"
  # print u.shape
  # print f.shape
  ##pre approx
  u_pre=GS_SOR(u,f,L,w,pre_step)
  r=residual(u_pre,f,L)
  if L>0:  
    r=restrict(r,L)
    ui=ini(L-1)
    for j in range(0,gamma):
      ux=mgm(ui,r,L-1,w,pre_step,post_step)
    u=u_pre+relax(ux,L-1)
  ##post approx
  u_post=GS_SOR(u,f,L,w,post_step)
  if L == 7 or L==8:
      r=residual(u_post,f,L)
      print "post {:d}".format(L)
      print (1./(2**L+1)**2*np.sum(np.abs(r)))
  return u_post

#Main 
#define parameters
start = time.time()
w=1.
pre_step=3
post_step=3
stepp=10000
gamma=1
L=8
N=2**L+2
ui=ini(L)
fxy=f(L)

err_plot=[]
n=[]

#Main multigrid iteration)
for t in range(1,10000):
  ux = mgm(ui,fxy,L,w,pre_step,post_step)
  r=residual(ux,fxy,L)
  err=(1./(N-2)**2*np.sum(np.abs(r)))
  print err
  err_plot.append(err)
  n.append(t)
  #print err
  if err<10**-8:
    break
  ui=ux

##Asymptatic fitting
# elapsed = (time.time() - start)
# err=np.array(err_plot)
# n=np.array(n)
# err=np.log(err)
# fit=np.polyfit(n[5:],err[5:],1)
# slope=fit[0]
# const=fit[1]
# xx=np.linspace(0,50,50)
# fitcurve=slope*xx+const
# print"w-circle, W={:f}".format(w)
# print"Time"
# print elapsed
# print"Slope"
# print slope


# plt.plot(err,"r+")
# plt.show()

#Plot
# converge, = plt.plot(n,err,"r+")
# fitting, = plt.plot(fitcurve,"b-")
# plt.legend([converge, fitting],['Experiment','Fitting'],loc=1)
# plt.axis([0,40,-30,5])
# plt.title("Convergence rate of V-circle Multi-grid algorithm (m1,m2)=(3,3) L=3")
# plt.xlabel("N circles")
# plt.ylabel("Log(Error)")
# plt.text(5,-20,"time={:.3f}\nLinear Fitting\ny={:.5f} x + {:.5f}".format(elapsed,slope,const))
# # plt.show()
# plt.savefig("V_converge_rate_L3.png")


##Plotting
# phi_m = np.linspace(0, 1, N)
# phi_p = np.linspace(0, 1, N)
# X,Y = np.meshgrid(phi_p, phi_m)
# exact=phi(L)

# fig = plt.figure(figsize=(14,6))

# # `ax` is a 3D-aware axis instance because of the projection='3d' keyword argument to add_subplot
# ax = fig.add_subplot(1, 2, 1, projection='3d')

# p = ax.plot_surface(X, Y, exact, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# cb = fig.colorbar(p, shrink=0.5)

# # surface_plot with color grading and color bar
# ax = fig.add_subplot(1, 2, 2, projection='3d')
# p = ax.plot_surface(X, Y, ux, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# cb = fig.colorbar(p, shrink=0.5)

# plt.savefig("Exact_Approx.png")

# print "multigrid method"
# print fit
# print n.shape



