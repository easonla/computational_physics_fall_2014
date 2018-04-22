##tau_int
from numpy import *
import matplotlib.pyplot as plt


def tau_int(u):
	n=len(u)
	sum=u[0]/2.
	for i in range(1,n):
		if i>5*sum:
			break
		else:
			sum=sum+u[i]
	return sum

#tau(Q,L,beta)
tau=zeros([3,5,6])
beta_all=array([0.1,0.2,0.3,0.4,0.5,1])
L_all=array([4,8,16,32,64])
a="Hot"
Nstep=120000

for i, L in enumerate(L_all):
	for j,Beta in enumerate(beta_all):
		print "loading file :{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt".format(a,L,Beta,Nstep) 
		data=loadtxt('acor/{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt'.format(a,L,Beta,Nstep))
		Ene=data[:,0]
		Mag=data[:,1]
		Mag_sq=data[:,2]
		tau[0,i,j]=tau_int(Ene)
		tau[1,i,j]=tau_int(Mag)
		tau[2,i,j]=tau_int(Mag_sq)


L_4, =  plt.semilogy(beta_all,tau[0,0,:],"ro-")
L_8, =  plt.semilogy(beta_all,tau[0,1,:],"bo-")
L_16, = plt.semilogy(beta_all,tau[0,2,:],"go-")
L_32, = plt.semilogy(beta_all,tau[0,3,:],"co-")
L_64, = plt.semilogy(beta_all,tau[0,4,:],"ko-")
plt.legend([L_4,L_8,L_16,L_32,L_64],['L=4','L=8','L=16','L=32','L=64'],loc=4)
plt.title(r"$\tau_{int}(\beta,L)$ of total energy")
plt.xlabel(r"$\beta$")
plt.ylabel(r'$\tau$')
plt.savefig("{:s}_Energy.png".format(a))
plt.show()

L_4, =  plt.semilogy(beta_all,tau[1,0,:],"ro-")
L_8, =  plt.semilogy(beta_all,tau[1,1,:],"bo-")
L_16, = plt.semilogy(beta_all,tau[1,2,:],"go-")
L_32, = plt.semilogy(beta_all,tau[1,3,:],"co-")
L_64, = plt.semilogy(beta_all,tau[1,4,:],"ko-")
plt.legend([L_4,L_8,L_16,L_32,L_64],['L=4','L=8','L=16','L=32','L=64'],loc=4)
plt.title(r"$\tau_{int}(\beta,L)$ of Magnetization")
plt.xlabel(r"$\beta$")
plt.ylabel(r'$\tau$')
plt.savefig("{:s}_Mag.png".format(a))
plt.show()

L_4, =  plt.semilogy(beta_all,tau[2,0,:],"ro-")
L_8, =  plt.semilogy(beta_all,tau[2,1,:],"bo-")
L_16, = plt.semilogy(beta_all,tau[2,2,:],"go-")
L_32, = plt.semilogy(beta_all,tau[2,3,:],"co-")
L_64, = plt.semilogy(beta_all,tau[2,4,:],"ko-")
plt.legend([L_4,L_8,L_16,L_32,L_64],['L=4','L=8','L=16','L=32','L=64'],loc=4)
plt.title(r"$\tau_{int}(\beta,L)$ of Magnetization")
plt.xlabel(r"$\beta$")
plt.ylabel(r'$\tau$')
plt.savefig("{:s}_Mag_sq.png".format(a))
plt.show()


