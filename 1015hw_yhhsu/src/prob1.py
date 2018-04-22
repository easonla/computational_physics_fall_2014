#computational physics hw5 
#numerical ODE, runge-kutta method
import math
import numpy as np
import matplotlib.pyplot as plt


def y_exact(t):
	return 1+eps*np.cos(t)
#Explicit Runge-kutta
def func(u0,u1,t):
	return u1
def dfunc(u0,u1,t):
	return -u0+const

step=100
#exact solution
ecc=0.1
const=1
phi=np.linspace(0,2.*np.pi,step)
u=1.+ecc*np.cos(phi)
r=1./u
x=r*np.cos(phi)
y=r*np.sin(phi)

#euler method
h=2.*np.pi/step
u_eu=np.zeros(step)
du_eu=np.zeros(step)
u_eu[0],du_eu[0]=1+ecc,0

for i in range(0,step-1):
	u_eu[i+1]=u_eu[i]+h*func(u_eu[i],du_eu[i],phi[i])
	du_eu[i+1]=du_eu[i]+h*dfunc(u_eu[i],du_eu[i],phi[i])

x_eu=1./u_eu*np.cos(phi)
y_eu=1./u_eu*np.sin(phi)



#RK2 method
def rk2(u0,u1,t):
	ak = dfunc(u0,u1+0.5*h*dfunc(u0,u1,t),t)
	u0 += h*u1
	u1 += h*ak
	return u0,u1

u_rk2=np.zeros(step)
du_rk2=np.zeros(step)
u_rk2[0],du_rk2[0]=1+ecc,0

for i in range(0,step-1):
	u_rk2[i+1],du_rk2[i+1]=rk2(u_rk2[i],du_rk2[i],phi[i])

x_rk2=1./u_rk2*np.cos(phi)
y_rk2=1./u_rk2*np.sin(phi)


#RK4 method
# def rk4(u0,u1,t):
# 	k1=dfunc(u0,u1,t)
# 	k2=dfunc(u0,u1+0.5*k1*h,t+0.5*h)
# 	k3=dfunc(u0,u1+0.5*k2*h,t+0.5*h)
# 	k4=dfunc(u0,u1+k3*h,t+h)
# 	u1+=h*(k1+2*k2+2*k3+k4)/6.

# 	k1=k2=k3=k4=u1
# 	u0+=h*(k1+2*k2+2*k3+k4)/6.

# 	return u0,u1
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

u_rk4=np.zeros(step)
du_rk4=np.zeros(step)
u_rk4[0],du_rk4[0]=1+ecc,0

for i in range(0,step-1):
	u_rk4[i+1],du_rk4[i+1]=rk4(u_rk4[i],du_rk4[i],phi[i])

x_rk4=1./u_rk4*np.cos(phi)
y_rk4=1./u_rk4*np.sin(phi)

#plotting
plt.plot(x_rk2,y_rk2,"r^")
plt.plot(x_rk4,y_rk4,"g^")
plt.plot(x_eu,y_eu,"b.")
plt.plot(x,y,"r-")
plt.show()