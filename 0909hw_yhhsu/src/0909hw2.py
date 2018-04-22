#9/9/2014 by yihsuan, hsu
#Computational Physics Class Homework 2
#given distance and lenght of an circular arc, find maximum height
import math

L = 1.0
S = 1.001
theta = math.sqrt(6*(1-L/S))
d =S/2/theta*(1-math.cos(theta))
print d