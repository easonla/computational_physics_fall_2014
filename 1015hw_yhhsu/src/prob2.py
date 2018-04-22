import math
import numpy as np
import matplotlib.pyplot as plt

def y_exact(t):
	return 1+ecc*np.cos(t)
def func(u0,u1,t):
	return u1
def dfunc(u0,u1,t):
	return -u0+const+3*lam*u0**2
def rk4(u0,u1,t):
	l1=u1
	k1=dfunc(u0,u1,t)
	l2=u1+0.5*h*k1
	k2=dfunc(u0+0.5*h*l1,u1+0.5*k1*h,t+0.5*h)
	l3=u1+0.5*h*k2
	k3=dfunc(u0+0.5*h*l2,u1+0.5*k2*h,t+0.5*h)
	l4=u1+h*k3
	k4=dfunc(u0+h*l3,u1+k3*h,t+h)

	u1+=h*(k1+2*k2+2*k3+k4)/6.
	u0+=h*(l1+2*l2+2*l3+l4)/6.

	return u0,u1
ecc_ar=np.array([0.1,0.3,0.5,0.9])
const=1
lam0=0.0001
name=(["epsilon=0.1","epsilon=0.3","epsilon=0.5","epsilon=0.9"])

for k in range (0,4):
	d_phi=[]
	lambdaa=[]
	ecc=ecc_ar[k]

	for p in range(1,100,10):
		lam=lam0*p
		u0=1.+ecc
		du0=0.
		step=10000
		h=2.5*np.pi/step

		phi=np.linspace(0,2.5*np.pi,step+1)
		u_exact=np.zeros(step+1)
		u_rk4=np.zeros(step+1)
		du_rk4=np.zeros(step+1)
		u_exact[0]=1
		u_rk4[0],du_rk4[0]=1+ecc,0

		for i in range(0,step):
			u_exact[i+1]=y_exact(phi[i+1])
			u_rk4[i+1],du_rk4[i+1]=rk4(u_rk4[i],du_rk4[i],phi[i])

		i=np.argmax(u_rk4[10:step+1])
		i+=10
		x0,x1,x2,x3,x4=u_rk4[i-1],u_rk4[i],u_rk4[i+1],u_rk4[i+2],u_rk4[i+3]
		y0,y1,y2=h*(i-1),h*(i),h*(i+1)
		d_phi.append(y1-2*np.pi)
		lambdaa.append(lam)

	lambdaaa=np.array(lambdaa)
	name[k],=plt.plot(lambdaa,d_phi)
	del lambdaa,d_phi

d_phi_ex=6*np.pi*lambdaaa
Linear,=plt.plot(lambdaaa,d_phi_ex,"k^-")
plt.legend([name[0],name[1],name[2],name[3],Linear],["$\epsilon=0.1$","$\epsilon=0.3$","$\epsilon=0.5$","$\epsilon=0.9$","Linear"],loc=4)
plt.xlabel('$\lambda$')
plt.ylabel('$\Delta \phi^{shift}$')
plt.savefig('problem2a.png')






