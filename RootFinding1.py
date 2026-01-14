import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
 

#Define Legendre Polynomials individually?



#Implementing Newton Raphson Method for root finding

#Recall we are working in the -1 to 1 domain 

#Also recall nth LP has n roots, need to find a way for the algorithm 
#to scan the entire domain 

'''
trial = legendre(2)
trialder = trial.deriv()

print(trial)

print(trialder)

print("see if I can eval ", trial(-0.57735))
'''

tol = 10**(-14)


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
			 
			RootArr.append(x0) #perhaps change this, its to find the next root
		
			print(x0)
	return RootArr

rootsP6 = NRalgorithm(6)

print("P0: ", NRalgorithm(0))
print("P1: ", NRalgorithm(1))
print("P2: ", NRalgorithm(2))
print("P3: ", NRalgorithm(3))

#print("Evaluating P1", f1(-0.57735))
#print("Second evaluation P1", f1(0.57735))
