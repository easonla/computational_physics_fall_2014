#9/9/2014 by yihsuan, hsu
#Computational Physics Class Homework 1 
#Numerical instability
import math
import numpy as np
from scipy.integrate import quad

# Exact Solution
# Reference 
# http://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html#general-integration-quad

def integrate(x,n):
	return x**n/(x+5)

# Approach solution
# Recursion formula I(n)=1/n-5*I(n-1)
def recursion(x,n):
	for i in range(1,n):
		x = 1.0/i - 5.0*x
	return x

#Main code

exact = []
exact_array = np.zeros((32,4,7,2))
for n in range(1,32):
	iter=quad(integrate,0,1,args=(n))
	#print iter[0]
	exact.append(iter[0])

single = np.float32(math.log(1.2))
double = np.float64(math.log(1.2))

#Output
print ("input\t"+ "initial_single = " +str(single)+ "\t" +"initial_double = " +str(double))
print ("n=\t" + "single value\t\t" + "double value\t\t" + "exact value\t\t" +"single error\t\t"+"double error\t\t")

for k in range(2,32):
	result_sin = recursion(single,k)
	result_dou = recursion(double,k)
	error_sin = exact[k-2] - result_sin
	error_dou = exact[k-2] - result_dou
	print str(k-1)+"\t\t"+str(result_sin)+"\t\t"+str(result_dou)+"\t\t"+str(exact[k-2])+"\t\t"+str(error_sin)+"\t\t"+str(error_dou)



f=open("recursion_table.txt","w")

f.write("input\t"+ "initial_single = " +str(single)+ "\t" +"initial_double = " +str(double)+"\n")
f.write("n=\t" + "single value\t\t" + "double value\t\t" + "exact value\t\t" +"single error\t\t"+"double error\t\t\n")

for k in range(2,32):
	result_sin = recursion(single,k)
	result_dou = recursion(double,k)
	error_sin = exact[k-2] - result_sin
	error_dou = exact[k-2] - result_dou
	f.write(str(k-1)+"\t\t"+str(result_sin)+"\t\t"+str(result_dou)+"\t\t"+str(exact[k-2])+"\t\t"+str(error_sin)+"\t\t"+str(error_dou)+"\n")


"{0:d} \t{1:.12e} \t{2:.12e}".format(k-1, result_sin, result_dou)

f.close()




