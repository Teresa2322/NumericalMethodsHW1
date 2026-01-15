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
	product = sp.prod(jarr)
	return sp.simplify(product)

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

plt.figure(1)
plt.plot(x_arr, P7(x_arr), label = 'P7', linestyle = '-') 

plt.plot(x_arr, phi_1f(x_arr), label = 'phi_1',  linestyle = 'dashed')
plt.plot(x_arr, phi_2f(x_arr), label = 'phi_2', linestyle = 'dashed')
plt.plot(x_arr, phi_3f(x_arr), label = 'phi_3', linestyle = 'dashed')
plt.plot(x_arr, phi_4f(x_arr), label = 'phi_4', linestyle = 'dashed')
plt.plot(x_arr, phi_5f(x_arr), label = 'phi_5', linestyle = 'dashed')
plt.plot(x_arr, phi_6f(x_arr), label = 'phi_6', linestyle = 'dashed')
plt.plot(x_arr, phi_7f(x_arr), label = 'phi_7', linestyle = 'dashed')

plt.plot(NRalgorithm(7), P7(NRalgorithm(7)), marker = 'o', color = 'b', ms = 5,  linestyle = 'None')

plt.xlabel("x")
plt.ylabel("y")
plt.legend()


#Computing the weights
#Will use simpson integration method, from class sample code
def simp( a, b, n, f):
	xk = np.linspace(a, b, 2*n+1)
	fk = f(xk)
	h = xk[1] - xk[0]
	return h*(fk[0] + fk[-1] + 4*fk[1:-1:2].sum() + 2*fk[2:-2:2].sum())/3.0

print("Trying out simpson", simp(-1,1,100,np.sin))

#function defining weights
def weights(N, k):
	phi = sp.lambdify(x, PolInd(N,k), modules='numpy')
	integral = simp(-1,1,5000,phi) #should probably reconsider how I def this 50000

	return integral 

print("Trying out weights", weights(7,1))

#Defining Quadrature Integration function

def QuadInt(N, Nsub, a, b, f):
	wk_arr  = []
	for i in range(1,N+1):
		wk_arr.append(weights(N, i)) #array of weights
	npwk_arr = np.array(wk_arr)
	xp_arr = np.array(NRalgorithm(N)) #array of sample points
	subint_arr = []
	for j in range(0,Nsub): #this will  go from 0 to Nsub-1 as needed
		
		#endpoints for a given subinterval 
		a_j = a + j*(b - a)/Nsub
		a_jp1 = a + (j+1)*(b - a)/Nsub
	
		#I believe we need to adjust the sample points and weights 
		#according to equations (2)-(4)
		xp_adj_arr = xp_arr*(a_jp1 - a_j)/2 + (a_jp1 + a_j)/2
		wk_adj_arr = npwk_arr*(a_jp1 - a_j)/2
		
		int_i = np.sum(wk_adj_arr*f(xp_adj_arr))
		subint_arr.append(int_i)
	return np.sum(subint_arr)

def f1(x):
	return np.exp(x)

f1_analytic = np.exp(1) - 1

def f2(x):
	return x**(1/3) + 1/(1+100*(x-5)**2)
# expected integral result f2 for given bounds: 15.7179202278389

f2_analytic = 15.7179202278389

def f3(x): 
	return (x**5)*np.abs(x)
# expected result for f3 for given bounds: − 431.4213979030774
f3_analytic = -431.4213979030776
print("testing quadrature function 1", QuadInt(7, 20, 0, 1, f1))
print("testing quadrature function 2", QuadInt(7, 20, 1, 10, f2)) 
print("Function 2 real value is:",  15.7179202278389)
print("testing quadrature function 3", QuadInt(7, 20, -np.pi, 4 - np.pi, f3))

#preliminary error analysis:

Err1_arr = []
Err2_arr = []
Err3_arr = []

for N_i in range(1,151):
	Evalf1_i = QuadInt(7, N_i, 0, 1, f1)
	Evalf2_i = QuadInt(7, N_i, 1, 10, f2)
	Evalf3_i = QuadInt(7, N_i, -np.pi, 4 - np.pi, f3)
	
	Err1_arr.append(np.abs(f1_analytic - Evalf1_i))
	Err2_arr.append(np.abs(f2_analytic - Evalf2_i))
	Err3_arr.append(np.abs(f3_analytic - Evalf3_i))
	

N_arr = np.linspace(1,150,150)

plt.figure(2)
plt.title("Function 1")
plt.loglog(N_arr*N_arr, Err1_arr)
plt.xlabel("N evaluations")
plt.ylabel("Abs Error")

plt.figure(3)
plt.title("Function 2")
plt.loglog(N_arr*N_arr, Err2_arr)
plt.xlabel("N evaluations")
plt.ylabel("Abs Error")

plt.figure(4)
plt.title("Function 3")
plt.loglog(N_arr*N_arr, Err3_arr)
plt.xlabel("N evaluations")
plt.ylabel("Abs Error")

plt.show()
