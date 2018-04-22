##Ising model, markov chain monte carlo and metropolis algorithm implement
from numpy import *
import matplotlib.pyplot as plt
import sys
import time
time_start=time.time()

starter=int(sys.argv[1])
Beta=float(sys.argv[2])
L=int(sys.argv[3])
if starter==1:
	a="Hot"
	print "Now is processing 'Hot' starter, Beta={:1.6f}, L={:d}".format(Beta,L)
else:
	a="Cold"
	print "Now is processing 'Cold' starter, Beta={:1.6f}, L={:d}".format(Beta,L)


def Hot(N):
	#Hot starter
	system=zeros((N,N),dtype=int)
	system=sign(2*random.rand(N,N)-1)
	return system

def Cold(N):
	#Cold starter
	system=zeros((N,N),dtype=int)
	system[:,:]=1
	return system

def CEnergy(system,i,j):
	N=len(system)
	return -1*system[i,j]*(system[(i+1)%N, j] + 
							system[i,(j+1)%N] + 
							system[(i-1)%N,j] + 
							system[i,(j-1)%N])

def TEnergy(system):
	#total energy
	N=len(system)
	totE=0
	for i in range(len(system)):
		for j in range(len(system)):
			totE+=-1*system[i,j]*(system[(i+1)%N, j] + 
									system[i,(j+1)%N] + 
									system[(i-1)%N,j] + 
									system[i,(j-1)%N])
	return totE/2.

def TEnergy_test(system):
	N=len(system)
	u=zeros([N+2,N+2])
	E=zeros([N,N])
	u[1:-1,-1]=system[:,0]
	u[1:-1,0]=system[:,-1]
	u[-1,1:-1]=system[0,:]
	u[0,1:-1]=system[-1,:]
	u[1:-1,1:-1]=system[:,:]
	E[:,:]= -1. * u[1:-1,1:-1] * (u[0:-2,1:-1] + u[2:,1:-1] + u[1:-1,0:-2] + u[1:-1,2:])	
	totE=sum(E)
	return totE/2.

def TMagnet(system):
	return sum(system)

def Montecarlo(system,N,T,Nstep):
	Iter=[]
	Ene=[]
	Mag=[]
	global count
	global step
	count=0		
	for step in range(0,Nstep):
		p = random.randint(0,N)
		q = random.randint(0,N)
		dE = -2.*CEnergy(system, p, q)

		if (dE <= 0.) or (dE > 0. and exp(-1./T*dE) > random.rand()):
			system[p,q] *= -1
			count+=1
		Iter.append(step)
		Ene.append(int(TEnergy(system)))
		Mag.append(int(TMagnet(system)))
		# if step%1000 ==1:
		# 	print step,TEnergy(system),TMagnet(system)	

	Iter=array(Iter)
	Ene=array(Ene)
	Mag=array(Mag)
	Mag_sq=Mag*Mag
	return Iter,Ene,Mag,Mag_sq

def Montecarlo_testing(system,N,T,Nstep):
	Iter=[]
	Ene=[]
	Mag=[]
	global count
	global step
	count=0		
	for step in range(0,Nstep):
		p = random.randint(0,N)
		q = random.randint(0,N)
		dE = -2.*CEnergy(system, p, q)

		if (dE <= 0.) or (dE > 0. and exp(-1./T*dE) > random.rand()):
			system[p,q] *= -1
			count+=1
		Iter.append(step)
		Ene.append(int(TEnergy_test(system)))
		Mag.append(int(TMagnet(system)))
		# if step%1000 ==1:
		# 	print step,TEnergy(system),TMagnet(system)	

	Iter=array(Iter)
	Ene=array(Ene)
	Mag=array(Mag)
	Mag_sq=Mag*Mag
	return Iter,Ene,Mag,Mag_sq


def time_plot(t,E,M,Ms):
	f, axarr = plt.subplots(3, sharex=True)
	axarr[0].plot(t,E,"-r")
	axarr[0].set_title('Energy $E(t)$')
	axarr[0].set_ylabel('$E$')
	axarr[1].plot(t,M,"-b")
	axarr[1].set_title('Magnetization  $M(t)$')
	axarr[1].set_ylabel('$M$')
	axarr[2].plot(t,Ms,"-g")
	axarr[2].set_title('Magnetization Square  $M^2(t)$')
	axarr[2].set_ylabel('$M^2$')
	axarr[2].set_xlabel("{:s} start. L={:d}, Beta={:1.6f} Iterations".format(a,Nside,Beta))
	plt.savefig("data/{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.png".format(a,Nside,Beta,Nstep))
	#plt.show()

def write_file(Iter,Ene,Mag,Mag_sq):
	f = open('data/{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt'.format(a,Nside,Beta,Nstep), 'w')
	f.write("##Iter	TotalE	TotalM\n")
	for i in range(0,Nstep):
		f.write("{:d} {:d} {:d} {:d}\n".format(Iter[i],Ene[i],Mag[i],Mag_sq[i]))
	f.close()	

if __name__ == '__main__':


	# Beta=1.
	Temp=1./Beta
	Nside=L
	Nstep=120000
	Niter=100000000
	if starter==1:
		system=Hot(Nside)
	else:
		system=Cold(Nside)

	# for i in range(0,20):
	# 	system=Hot(Nside)
	# 	print TEnergy(system)
	# 	print TEnergy_test(system)

	system1=system
	t=time.time()
	Iter,Ene,Mag,Mag_sq=Montecarlo_testing(system,Nside,Temp,Nstep)
	accep_rate=float(count)/float(step)
	write_file(Iter,Ene,Mag,Mag_sq)
	time_plot(Iter,Ene,Mag,Mag_sq)
	end=time.time()
	eclipse=end - time_start
	print "Acceptence rate = {:1.3f}, Time eclipse={:f}".format(accep_rate,eclipse)
	print "="*70




