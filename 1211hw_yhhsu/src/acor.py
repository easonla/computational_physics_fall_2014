#autocorrelation time, autocorrelation function
from numpy import *
import matplotlib.pyplot as plt
import sys
import time
time_start=time.time()
def write_file(autocorrelation):
	f = open('acor/{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt'.format(a,L,Beta,Nstep), 'w')
	f.write("##autocorrelation function\n")
	f.write("##{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt\n".format(a,L,Beta,Nstep))
	j=len(autocorrelation[:,1])
	for i in range(0,j):
		f.write("{:f} {:f} {:f}\n".format(autocorrelation[i,0],autocorrelation[i,1],autocorrelation[i,2]))
	f.close()

def plot(autocorrelation):
	f, axarr = plt.subplots(3, sharex=True)
	axarr[0].plot(autocorrelation[:,0],"-r")
	axarr[0].set_ylabel('$\\rho_E(t)$')
	axarr[1].plot(autocorrelation[:,1],"-b")
	axarr[1].set_ylabel('$\\rho_M(t)$')
	axarr[2].plot(autocorrelation[:,2],"-g")
	axarr[2].set_ylabel('$\\rho_{M^2}(t)$')
	axarr[2].set_xlabel("{:s} start. L={:d}, Beta={:1.6f} steps".format(a,L,Beta))
	plt.savefig("acor/{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.png".format(a,L,Beta,Nstep))

starter=int(sys.argv[1])
Beta=float(sys.argv[2])
L=int(sys.argv[3])
Nstep=120000
if starter==1:
	a="Hot"
	#print "'Hot' starter, Beta={:1.6f}, L={:d}".format(Beta,L)
else:
	a="Cold"
	#print "'Cold' starter, Beta={:1.6f}, L={:d}".format(Beta,L)

print "loading file :{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt".format(a,L,Beta,Nstep) 
data=loadtxt('data/{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt'.format(a,L,Beta,Nstep), dtype=int)

pre=0
Ene=data[pre:,1]
Mag=data[pre:,2]
Mag_sq=data[pre:,3]
nstep=len(Ene)

time_test=time.time()
k=120000
autocorrelation=zeros((k,3))
for i, Q in enumerate([Ene,Mag,Mag_sq]):
	Mean=sum(Q)/nstep
	Mean_sq=sum(Q**2)/nstep
	for t in range(0,k):
		acorr=0.
		acorr=dot(Q[:nstep-t].T,Q[t:nstep])
		Mean_corr=acorr/float(nstep-t)
		autocorrelation[t,i]=(Mean_corr-Mean**2)/(Mean_sq-Mean**2)
time_test_end=time.time()
time_sm=time_test_end-time_test
plot(autocorrelation)
write_file(autocorrelation)
print "="*20+"time=",time_sm,"="*20











