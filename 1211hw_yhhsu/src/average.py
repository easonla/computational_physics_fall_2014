from numpy import *
import matplotlib.pyplot as plt

a="Hot"
L=16
Nstep=120000
cutoff=20000
beta_all=array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
Ene_av=zeros(10)
Mag_av=zeros(10)
Mag_sq_av=zeros(10)

for j,Beta in enumerate(beta_all):
	print "loading file :{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt".format(a,L,Beta,Nstep) 
	data=loadtxt('data/{:s}_L_{:d}_Beta_{:1.6f}_iter_{:d}.txt'.format(a,L,Beta,Nstep))
	Ene=data[:,1]
	Mag=data[:,2]
	Mag_sq=data[:,3]
	Ene_av[j]=sum(Ene[cutoff:])/float(Nstep-cutoff)
	Mag_av[j]=sum(Mag[cutoff:])/float(Nstep-cutoff)
	Mag_sq_av[j]=sum(Mag_sq[cutoff:])/float(Nstep-cutoff)


f, axarr = plt.subplots(3, sharex=True)
axarr[0].plot(beta_all,Ene_av,"-or")
axarr[0].set_ylabel(r'$<E(t)>$')
axarr[1].plot(beta_all,Mag_av,"-ob")
axarr[1].set_ylabel(r'$<M(t)>$')
axarr[2].plot(beta_all,Mag_sq_av,"-og")
axarr[2].set_ylabel(r'$<M(t)^2>$')
axarr[2].set_xlabel("{:s} start. L={:d} 120000 steps ,20000 cutoff steps".format(a,L))
plt.savefig("average_plot_L_16_Hot.png")
plt.show()