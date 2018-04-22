import math

# #python exercise
# #list is an object that can store list of different variation
# a = 4.4234
# b = 10
# c = a+b 
# list = [a,b,c,"abc","i miss taiwan"]
# #print list
# list
# #print list4
# list[3]



#define a new function
def myexp(x, n):
	sum  = 0
	term = 1.0



	for i in range(n+1):
		
		sum = sum + term
		term = term * x/(i+1)

	return sum


x1 = 20
y1 = myexp(x1,10)
print(y1)

y_exact= math.exp(x1)
print(y_exact)

print("\nour function"= str(y2))
print("exact function"= str(y_exact)+"\n")