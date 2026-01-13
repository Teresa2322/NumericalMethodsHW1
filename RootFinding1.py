import numpy as np

#Define Legendre Polynomials individually?



#Implementing Newton Raphson Method for root finding

#Recall we are working in the -1 to 1 domain 

#Also recall nth LP has n roots, need to find a way for the algorithm 
#to scan the entire domain 


tol = 10**(-14)

def P1(x):
	return (3*x**2 - 1)/2

def grad(
def NRalgorithm( x, f(x), derf(x) ):
	x0 = -1 
	RootArr = []
	while abs(f(x0)) > tol and x0 < 1:
		x0 = x0 - (f(x0)/derf(x1))
	else: 
		RootArr.append(x0)
		x0 = x0 + 0.00001 #perhaps change this, its to find the next root

	return RootArr

