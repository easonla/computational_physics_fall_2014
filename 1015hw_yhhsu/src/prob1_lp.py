#computational physics hw5 
#numerical ODE, runge-kutta method
import math
import numpy as np
import matplotlib.pyplot as plt

ecc=0.9
const=1
u0=1+ecc
du0=0
ss=[]
er_eu=[]
der_eu=[]
er_rk2=[]
der_rk2=[]
er_rk4=[]
der_rk4=[]

def y_exact(t):
	return 1+ecc*np.cos(t)
def func(u0,u1,t):
	return u1
def dfunc(u0,u1,t):
	return -u0+const
def euler(u0,u1,t):
	u0+=h*func(u0,u1,t)
	u1+=h*dfunc(u0,u1,t)
	return u0,u1
def rk2(u0,u1,t):
	l1=u1
	k1=dfunc(u0,u1,t)
	l2=u1+0.5*h*k1
	k2=dfunc(u0+0.5*h*l1,u1+0.5*k1*h,t+0.5*h)
	u0 += h*l2
	u1 += h*k2
	return u0,u1
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

for step in range(50,1000,10):
	h=2.*np.pi/step
	ss.append(h)
	phi=np.linspace(0,2.*np.pi,step+1)

	u_eu=np.zeros(step+1)
	du_eu=np.zeros(step+1)
	u_eu[0],du_eu[0]=1+ecc,0

	u_rk2=np.zeros(step+1)
	du_rk2=np.zeros(step+1)
	u_rk2[0],du_rk2[0]=1+ecc,0

	u_rk4=np.zeros(step+1)
	du_rk4=np.zeros(step+1)
	u_rk4[0],du_rk4[0]=1+ecc,0

	for i in range(0,step):
		u_eu[i+1],du_eu[i+1]=euler(u_eu[i],du_eu[i],phi[i])
		u_rk2[i+1],du_rk2[i+1]=rk2(u_rk2[i],du_rk2[i],phi[i])
		u_rk4[i+1],du_rk4[i+1]=rk4(u_rk4[i],du_rk4[i],phi[i])

	er_eu.append(math.fabs(u_eu[step])-u0)
	der_eu.append(math.fabs(du_eu[step]))

	er_rk2.append(math.fabs(math.fabs(u_rk2[step])-u0))
	der_rk2.append(math.fabs(du_rk2[step]))
	
	er_rk4.append(math.fabs(math.fabs(u_rk4[step])-u0))
	der_rk4.append(math.fabs(du_rk4[step]))

#plotting
# u=1.+ecc*np.cos(phi)
# r=1./u
# x=r*np.cos(phi)
# y=r*np.sin(phi)
# x_eu=1./u_eu*np.cos(phi)
# y_eu=1./u_eu*np.sin(phi)
# x_rk2=1./u_rk2*np.cos(phi)
# y_rk2=1./u_rk2*np.sin(phi)
# x_rk4=1./u_rk4*np.cos(phi)
# y_rk4=1./u_rk4*np.sin(phi)

# plt.plot(x_rk2,y_rk2,"b-")
# plt.plot(x_rk4,y_rk4,"g^")
# plt.plot(x_eu,y_eu,"ro")
# plt.plot(x,y,"r-")
# plt.show()
eu_co=np.zeros(step+1)
rk2_co=np.zeros(step+1)
rk4_co=np.zeros(step+1)

eu_co=(1/np.array(ss))**-1
rk2_co=(1/np.array(ss))**-2
rk4_co=(1/np.array(ss))**-4

deu_co=np.zeros(step+1)
drk2_co=np.zeros(step+1)
drk4_co=np.zeros(step+1)

deu_co=(1/np.array(ss))**-1
drk2_co=(1/np.array(ss))**-2
drk4_co=(1/np.array(ss))**-4

plt.subplot(211)
Euler_Mtod,=plt.loglog(ss,er_eu,"k.-")
RK2_Mtod,=plt.loglog(ss,er_rk2,"b.-")
RK4_Mtod,=plt.loglog(ss,er_rk4,"r.-")
Euler_co,=plt.loglog(ss,eu_co,"k--")
RK2_co,=plt.loglog(ss,rk2_co,"b--")
RK4_co,=plt.loglog(ss,rk4_co,"r--")
plt.legend([Euler_Mtod,RK2_Mtod,RK4_Mtod],['Euler_Method','RK2_Method','RK4_Method'],loc=4)
plt.title("error in $u$")
plt.xlabel('step size h')
plt.ylabel('$\epsilon_N$')

plt.subplot(212)
dEuler_Mtod,=plt.loglog(ss,der_eu,"k.-")
dRK2_Mtod,=plt.loglog(ss,der_rk2,"b.-")
dRK4_Mtod,=plt.loglog(ss,der_rk4,"r.-")
dEuler_co,=plt.loglog(ss,deu_co,"k--")
dRK2_co,=plt.loglog(ss,drk2_co,"b--")
dRK4_co,=plt.loglog(ss,drk4_co,"r--")
plt.legend([dEuler_Mtod,dRK2_Mtod,dRK4_Mtod],['Euler_Method','RK2_Method','RK4_Method'],loc=4)
plt.title("error in $du/d(\phi)$")
plt.xlabel('step size h')
plt.ylabel('$\epsilon_N$')

plt.suptitle('Eccentricity = 0.9',fontsize=18)
plt.savefig('error_ecc09.png')
plt.show()










