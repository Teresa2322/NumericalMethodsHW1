import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
import sympy as sp 

#Define Legendre Polynomials individually?



#Implementing Newton Raphson Method for root finding

#Recall we are working in the -1 to 1 domain 

#Also recall nth LP has n roots, need to find a way for the algorithm 
#to scan the entire domain 

tol = 10**(-15)
x = sp.symbols('x')
x_arr = np.linspace(-1,1,1000)

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

#Finding the indicator polynomials

def PolInd(N, k):
	xpi = NRalgorithm(N)
	jarr = []
	for i in range(1,N+1):
		if i != k: #end point in range is otherwise excluded
			j = (x - xpi[i-1])/(xpi[k-1] - xpi[i-1]) #accounting for index starts at 0
			jarr.append(j)
	product = np.prod(jarr)
	return product

print("Trial Polynomial", PolInd(3,1))
print("Trial Evaluated Polynomial", PolInd(3,1).subs(x, 0.5))

#N = 7 Legendre Polynomial
P7 = legendre(7)

#indicator polynomials for N = 7
phi_1f = sp.lambdify(x, PolInd(7,1), modules='numpy')
phi_2f = sp.lambdify(x, PolInd(7,2), modules='numpy')
phi_3f = sp.lambdify(x, PolInd(7,3), modules='numpy')
phi_4f = sp.lambdify(x, PolInd(7,4), modules='numpy')
phi_5f = sp.lambdify(x, PolInd(7,5), modules='numpy')
phi_6f = sp.lambdify(x, PolInd(7,6), modules='numpy')
phi_7f = sp.lambdify(x, PolInd(7,7), modules='numpy')

#plot


plt.plot(x_arr, P7(x_arr), label = 'P7') 
plt.plot(x_arr, phi_1f(x_arr), label = 'phi_1')
plt.plot(x_arr, phi_2f(x_arr), label = 'phi_2')
plt.plot(x_arr, phi_3f(x_arr), label = 'phi_3')
plt.plot(x_arr, phi_4f(x_arr), label = 'phi_4')
plt.plot(x_arr, phi_5f(x_arr), label = 'phi_5')
plt.plot(x_arr, phi_6f(x_arr), label = 'phi_6')
plt.plot(x_arr, phi_7f(x_arr), label = 'phi_7')
plt.legend()
plt.show()

