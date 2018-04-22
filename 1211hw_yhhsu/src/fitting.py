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
Chi=zeros(7)

for i, L in enumerate(L_all):
	print "loading file :{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt".format(a,L,Beta,Nstep) 
	data=loadtxt('data/{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt'.format(a,L,Beta,Nstep), dtype=int)
	N=len(data[pre:,0])
	Mag_sq=data[pre:,3]
	Mean_sq=sum(Mag_sq)/N
	Chi[i]=Mean_sq/float(L**2)


Chi_log=zeros(7)
Chi_log=log(Chi)
L_all_log=zeros(7)
L_all_log=log(L_all)

fit=polyfit(L_all_log[:4],Chi_log[:4],1)
slope=fit[0]
const=fit[1]
xx=linspace(4,256)
fitcurve=exp(const)*xx**slope

#Plot
converge, = plt.loglog(L_all,Chi,"ro")
fitting, = plt.loglog(xx,fitcurve,"b-")
plt.legend([converge, fitting],['Experiment','Fitting'],loc=4)
plt.title(r"susceptibilitiy $\chi$ versus $L$ ,{:s} start".format(a))
plt.xlabel(r"$L$")
plt.ylabel(r"$\chi$")
plt.text(10,1000,r"$\chi$ ~ $L^\alpha$, $\alpha$ = {:1.2f}".format(slope))
plt.savefig("{:s}_susceptability.png".format(a))
plt.show()
