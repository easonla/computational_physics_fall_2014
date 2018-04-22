#9/30 Hw3 yihsuan, hsu
#searching for bifurcation point
import math
import random
import sys
import numpy as np
import matplotlib.pyplot as plt

##define querty function
"""
def func(a,x):
	return a*x*(1.-x)
def dfunc(a,x):
	return a*(1.-2.*x)
"""
def func(a,x):
	return a*math.sin(math.pi*x)
def dfunc(a,x):
	return a*math.pi*math.cos(math.pi*x)


##defin m-degree iterated function and its derivative
def m_dfunc(a,m,x):
	dx=1
	for i in range(0,m+1):
		m_dx = dfunc(a,x)
		x = func(a,x)
		dx=m_dx*dx
	return dx
def m_func(a,m,x):
	for j in range (0,2**m):
		x = func(a,x)
	return x
##define logistic
# parameter a, circle m, initial value x
def logistic(a,m,x):
	for i in range (1,10000):
		xi = x
		for j in range (0,2**m):
			x = func(a,x)
	return x
def Newton(a,m,x):
	for i in range(0,10002):
		xi=x-(m_func(a,m,x)-x)/(m_dfunc(a,m,x)-1.)
		if math.fabs(xi-x)<10**-10:
			#print "got it, i={}".format(i)
			break
		x=xi
		if i>10000:
			xi=0
	return xi

# print logistic(3.1,0,0.4)
# print m_dfunc(3.,0,0.4),m_func(3.,0,0.4)

#sys.argv user input
#Main
"""
n	cycle (2^n)	r_(2^n)	
1	2	3	
2	4	3.449490	
3	8	3.544090	
4	16	3.564407	
5	32	3.568750	
6	64	3.56969	
7	128	3.56989	
8	256	3.569934	
9	512	3.569943	
10	1024	3.5699451	
11	2048	3.569945557	
infty	accumulation point	
"""

method=int(sys.argv[1])
m=int(sys.argv[2])
a=float(sys.argv[3])
h=float(sys.argv[4])
##hand adjustment
x0=random.random()

if method==1:
	#for k in range(10000):
	for i in range(1000000):
		a=a+i*10**(-h)
		xn=logistic(a,m,x0)
		x_ex=Newton(a,m,xn)
		lam=m_dfunc(a,m,x_ex)
		if x_ex==0:
			x0=random.random()
			print "Newton Fail change initial value = {}".format(x0)
		elif lam<-1:
			print "Success :lam={:.12f} a={:.13f} step={:4d} accuracy={:g}".format(lam,a,i,10**(-h))
			break
		elif i%10==1:
			print "Fail    :lam={:.13g} a={:.13f} step={:4d}".format(lam,a,i)
		else:
			if x_ex==0:
				x0=random.random()

##bisection
#{80.31831, 0.719962, 0.833266, 0.858609, 0.864084, 0.865259}

if method==2:
#	a_l,a_up=a,a+1000.*10**(-h)
#	for k in range(100):
#	for f in range(0,50):
	acc=10**-8
	a_l,a_up=0.86,0.87
	x0=random.random()
	lam_test=[]
	a_test=[]
	print "2^{}-period searching ".format(m)
	print "initial x0={:.5f}, initial bracket=[{:.5f},{:.5f}]".format(x0,a_l,a_up)
	for i in range(300):
		#x0=random.random()
		#lam_test.append(lam)
		#a_test.append(a)
		a=0.5*(a_l+a_up)
		xn=logistic(a,m,x0)
		x_ex=Newton(a,m,xn)
		lam=m_dfunc(a,m,x_ex)


		if math.fabs(1.+lam) < acc:
			print "Success :lam={:.12f} a={:.13f} step={:4d} accuracy={:g}".format(lam,a,i,acc)
			break
		elif i%10==1:
			print "Fail    :lam={:.13g} a={:.13f} step={:4d}".format(lam,a,i)
		else:
			if x_ex==0:
				x0=random.random()
				print "Newton fail, change initial value = {:.5f}".format(x0)
			elif lam<-10**6:
				x0=random.random()
				print "Lambda fail, change initial value = {:.5f}".format(x0)
			# elif lam-lam_test<10**-5:## reset if run out of range
			# 	a_l,a_up,x0=3.56,3.567,random.random()
			elif lam<-1.:
				a_up=a
			else:
				a_l=a

"""
	plt.plot(a_test,lam_test,'ro')
	plt.axis([3.44814,3.44816,-1-10**-5,-1+10**-5])
	plt.grid(True)
	plt.show()
"""

"""
	from itertools import permutations
	from random import sample
	Nlines = 200
	color_lvl = 8
	rgb = np.array(list(permutations(range(0,256,color_lvl),3)))/255.0
	for k in range(0,200):
		print rgb[k]
	colors = sample(rgb,Nlines)
	for k in range(0,200):
		print colors[k]

	for k in range(0,i):
		plt.plot(a_test[k],lam_test[k],'^',color=colors[k])
	plt.show()
"""


