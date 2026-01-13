import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
 

#Define Legendre Polynomials individually?



#Implementing Newton Raphson Method for root finding

#Recall we are working in the -1 to 1 domain 

#Also recall nth LP has n roots, need to find a way for the algorithm 
#to scan the entire domain 


trial = legendre(2)
trialder = trial.deriv()

print(trial)

print(trialder)

print("see if I can eval ", trial(-0.57735))

tol = 10**(-14)


def NRalgorithm( x, n):
	
	x0 = -1 
	RootArr = []
	while abs(f(x0)) > tol and x0 < 1:
		x0 = x0 - (f(x0)/derf(x1))
	else: 
		RootArr.append(x0)
		x0 = x0 + 0.00001 #perhaps change this, its to find the next root

	return RootArr


#print("Evaluating P1", f1(-0.57735))
#print("Second evaluation P1", f1(0.57735))
