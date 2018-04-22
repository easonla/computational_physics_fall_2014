#9/18/2014 by yihsuan hsu
#computational physics hw 2
#one dimension root finding through bisection, regular falsi, secant and Newton method
import math
import matplotlib.pyplot as plt
import numpy as np

#bracket, iteration, function
up,low=3.0,3.5
#up,low=np.float64(3.0),np.float64(3.5)
iteration=100000000
exact_root=math.pi
xacc=1e-14

def func(x):
	#y=(x-1.0)*(x-10.0)*(x-100.0)
	# y= x-100
	# y = math.log(x,2)
	y = math.sin(x)
	# y=math.exp(x)
	# y=1/(x-1)
	# y=x**11
	# y=x*math.e**(-1/x**2)
	return y
def dif(x):
	return (func(x+1e-14)-func(x))/1e-14

def stop_func(x):
	return test

def bisection(a,b,n):
	itr=[]
	root=[]
	eps=[]
	xn,xn_1=0,0

	for i in range (0,n):
		mid = (a + b)/2.0
		if math.fabs(mid-exact_root)>xacc:
			root.append(mid)
			itr.append(math.fabs(func(mid)))
			if func(a)*func(mid)>0:
				a=mid
			else:
				b=mid

			xn=mid
			eps.append(math.fabs(xn-xn_1))
			xn_1=mid
		else:
			break

	return i,itr,root,eps

def regula_falsi(a,b,n):
	itr=[]
	root=[]
	eps=[]
	xn,xn_1=0,0
	for i in range (0,n):
		w = (func(b)*a-func(a)*b)/(func(b)-func(a))
		if math.fabs(w-exact_root)>xacc:
			root.append(w)
			itr.append(math.fabs(func(w)))
			if func(a)*func(w)<=0:
				eps.append(math.fabs(b-w))
				b=w
			else:
				eps.append(math.fabs(a-w))
				a=w

			xn=w
			eps.append(math.fabs(xn-xn_1))
			xn_1=w
		else:
			break
	return i,itr,root,eps

def secant(a,b,n):
	itr=[]
	root=[]
	eps=[]
	xn,xn_1=0,0
	for i in range (0,n):
		w = a - (func(a)*(a-b))/(func(a)-func(b))
		if math.fabs(w-exact_root)>xacc:
			root.append(w)
			itr.append(math.fabs(func(w)))

			if func(a)*func(w)<=0:
				b=w
			else:
				a=w	

			xn=w
			eps.append(math.fabs(xn-xn_1))
			xn_1=w
		else:
			break

	return i,itr,root,eps

def newton(a,b,n):
	itr=[]
	root=[]
	eps=[]
	for i in range (0,n):
		w = a - func(a)/dif(a)
		if math.fabs(w-exact_root)>xacc:
			itr.append(func(w))
			root.append(w)
			eps.append(math.fabs(a-w))
			a = w
		else:
			break
	return i,itr,root,eps

#main code

if func(up)*func(low)<0:
	n_bi,bi_func,bi_root,bi_eps = bisection(up,low,iteration)
	n_rf,rf_func,rf_root,rf_eps = regula_falsi(up,low,iteration)
	n_se,se_func,se_root,se_eps = secant(up,low,iteration)
	n_nt,nt_func,nt_root,nt_eps = newton(up,low,iteration)
else:
	print "bracket is invalid"
	sys.exit()

# up,low=3.1,3.2
# if func(up)*func(low)<0:
# 	bi_func_1,bi_root_1,bi_eps_1 = bisection(up,low,iteration)
# 	rf_func_1,rf_root_1,rf_eps_1 = regula_falsi(up,low,iteration)
# 	se_func_1,se_root_1,se_eps_1 = secant(up,low,iteration)
# 	nt_func_1,nt_root_1,nt_eps_1 = newton(up,low,iteration)
# else:
# 	print "bracket is invalid"
# 	sys.exit()


#output
bi_err=[]
rf_err=[]
se_err=[]
nt_err=[]

#print "\n bisection iteration \n"
#print "iteration\tf(Xn)\troot\tepsilon\terror"
for i in range(0,n_bi):
	bi_err.append(math.fabs(exact_root-bi_root[i]))
#	print i,bi_func[i],bi_root[i],bi_eps[i],bi_err[i]

#print "\n ragular falsi iteration \n"
#print "iteration\tf(Xn)\troot\tepsilon\terror"
for i in range(0,n_rf):
	rf_err.append(math.fabs(exact_root-rf_root[i]))
#	print i,rf_func[i],rf_root[i],rf_eps[i],rf_err[i]

#print "\n secant method \n"
#print "iteration\tf(Xn)\troot\tepsilon\terror"
for i in range(0,n_se):
	se_err.append(math.fabs(exact_root-se_root[i]))
#	print i,se_func[i],se_root[i],se_eps[i],rf_err[i]

#print "\n newton method \n"
#print "iteration\tf(Xn)\troot\tepsilon\terror"
for i in range(0,n_nt):
	nt_err.append(math.fabs(exact_root-nt_root[i]))
	#print i,nt_func[i],nt_root[i],nt_eps[i],nt_err[i]

#make plots

plt.figure()
plt.subplot(2,2,1)
plt.axis([0,n_bi,1e-16,1])
plt.plot(bi_func,'ro')
plt.semilogy()
plt.xlabel("bisection")
plt.ylabel("f(Xn)")

plt.subplot(2,2,2)
plt.axis([0,n_rf,1e-16,1])
plt.semilogy()
plt.plot(rf_func,'ro')
plt.xlabel("ragular falsi")

plt.subplot(2,2,3)
plt.plot(se_func,'ro')
plt.axis([0,n_se,1e-16,1])
plt.semilogy()
plt.xlabel("secant")

plt.subplot(2,2,4)
plt.plot(nt_func,'ro')
plt.axis([0,n_nt,1e-16,1])
plt.semilogy()
plt.xlabel("newton")

plt.suptitle("numerical root finding sin(x)^3 in [3.1,3.2]")
plt.savefig("func_sin_tri_1.png")

plt.show()


plt.figure()
plt.subplot(2,2,1)
plt.axis([0,n_bi,1e-16,1])
plt.plot(bi_eps,'ro')
plt.semilogy()
plt.xlabel("bisection")
plt.ylabel("epsilon = |x_n - x_n-1|")

plt.subplot(2,2,2)
plt.axis([0,n_rf,1e-16,1])
plt.semilogy()
plt.plot(rf_eps,'ro')
plt.xlabel("ragular falsi")

plt.subplot(2,2,3)
plt.plot(se_eps,'ro')
plt.axis([0,n_se,1e-16,1])
plt.semilogy()
plt.xlabel("secant")

plt.subplot(2,2,4)
plt.plot(nt_eps,'ro')
plt.axis([0,n_nt,1e-16,1])
plt.semilogy()
plt.xlabel("newton")

plt.suptitle("numerical root finding sin(x) [3.0,3.5] in red by autostop algorithm")
plt.savefig("eps_sin_auto.png")

plt.show()

plt.figure()
plt.subplot(2,2,1)
plt.axis([0,n_bi,1e-16,1])
plt.plot(bi_err,'ro')
plt.semilogy()
plt.xlabel("bisection")
plt.ylabel("err = |x_n - x_exact|")

plt.subplot(2,2,2)
plt.axis([0,n_rf,1e-16,1])
plt.semilogy()
plt.plot(rf_err,'ro')
plt.xlabel("ragular falsi")

plt.subplot(2,2,3)
plt.plot(se_err,'ro')
plt.axis([0,n_se,1e-16,1])
plt.semilogy()
plt.xlabel("secant")

plt.subplot(2,2,4)
plt.plot(nt_err,'ro')
plt.axis([0,n_nt,1e-16,1])
plt.semilogy()
plt.xlabel("newton")

plt.suptitle("numerical root finding sin(x) in [3.0,3.5]")
plt.savefig("err_sin_auto.png")

plt.show()
