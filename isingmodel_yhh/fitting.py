###least square fitting
from numpy import *
import scipy.optimize as optimization
import matplotlib.pyplot as plt

"""
def func(x, a, b, c):
    return a + b*x + c*x*x

# Generate artificial data = straight line with a=0 and b=1
# plus some noise.
xdata = numpy.array([0.0,1.0,2.0,3.0,4.0,5.0])
ydata = numpy.array([0.1,0.9,2.2,2.8,3.9,5.1])
# Initial guess.
x0    = numpy.array([0.0, 0.0, 0.0])

sigma = numpy.array([1.0,1.0,1.0,1.0,1.0,1.0])

print optimization.curve_fit(func, xdata, ydata, x0, sigma)
"""

starter=1
Beta=0.440687
L_all=array([4,8,16,32,64,128,256])
Nstep=120000
pre=20000
a="Hot"
Kai=zeros(7)

for i, L in enumerate(L_all):
	print "loading file :{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt".format(a,L,Beta,Nstep) 
	data=loadtxt('data/{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt'.format(a,L,Beta,Nstep), dtype=int)
	N=len(data[pre:,0])
	Mag_sq=data[:,3]
	Mean_sq=sum(Mag_sq)/N
	Kai[i]=Mean_sq/L**2

plt.loglog(L_all,Kai,"o")
plt.show()