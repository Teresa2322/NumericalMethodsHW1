import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
 

#Define Legendre Polynomials individually?



#Implementing Newton Raphson Method for root finding

#Recall we are working in the -1 to 1 domain 

#Also recall nth LP has n roots, need to find a way for the algorithm 
#to scan the entire domain 

tol = 10**(-15)


def NRalgorithm(n):

	func = legendre(n)
	derf = func.deriv()

	x0 = -1 
	RootArr = []	

	while x0 < 1:
		x0 += 0.01
		if func(x0)*func(x0 + 0.01) < 0: #using this to id a region with a root  
			while abs(func(x0)) > tol:
				x0 = x0 - (func(x0)/derf(x0))			 
			RootArr.append(x0) 
	return RootArr

#collecting sample points (in arrays, by degree of legendre polynomial)
xp0 = NRalgorithm(0)
xp1 = NRalgorithm(1)
xp2 = NRalgorithm(2)
xp3 = NRalgorithm(3)
xp4 = NRalgorithm(4)
xp5 = NRalgorithm(5)
xp6 = NRalgorithm(6)
xp7 = NRalgorithm(7)

#Sanity check: checking number of roots:
print("length checks:", len(xp0), len(xp1), len(xp2), len(xp3), len(xp4), len(xp5), len(xp6), len(xp7))


