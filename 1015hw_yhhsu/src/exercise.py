#Computational physics, Runge-kutta method, adaptive step-size control
#multi-step methods, predictor correction, stability, stiff equation Exercise.
import math
import numpy as np
import matplotlib.pyplot as plt

#simple euler methods
#y'=y(1-y) , y(0)=1/10
#y_exact = 1 / (9*exp(-t)+1)

def euler1(y):
	return y*(1-y)
def y_exact(t):
	return 1/(9*math.exp(-t)+1)

#Explicit Runge-kutta
def func(y,t):
	return y*(1-y)
	u0
def rk2(y,t):
	ak = func(y+0.5*h*func(y,t),t+0.5*h)
	return h*ak

def rk4(y,t):
	k1 = func(y,t)
	k2 = func(y+0.5*h*k1,t+0.5*h)
	k3 = func(y+0.5*h*k2,t+0.5*h)
	k4 = func(y+h*k3,t+h)
	return 1./6.*(k1+2*k2+2*k3+k4)*h


h = 2
y0=0.1
num = 15
y_ex = np.zeros(num+1)
yn_eu = np.zeros(num+1)
yn_rk2 = np.zeros(num+1)
yn_rk4 = np.zeros(num+1)
t = np.linspace(0,h*num,num+1)

#initial value 
y_ex[0] = y_exact(t[0])
yn_eu[0]=y0
yn_rk2[0]=y0
yn_rk4[0]=y0

for i in range(0,num):
	yn_eu[i+1]=yn_eu[i]+h*func(yn_eu[i],t[i])
	yn_rk2[i+1]=yn_rk2[i]+rk2(yn_rk2[i],t[i])
	yn_rk4[i+1]=yn_rk4[i]+rk4(yn_rk4[i],t[i])
	y_ex[i+1] = y_exact(t[i+1])

plt.axis([0,30,0,2])
plt.plot(t,yn_eu,"b-")
plt.plot(t,y_ex,"r-")
plt.plot(t,yn_rk2,"g-")
plt.plot(t,yn_rk4,"y-")
plt.show()