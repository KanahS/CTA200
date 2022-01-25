# Program to calculate binary evolution from constant time lag formalism

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from scipy.integrate import odeint
from scipy.integrate import quad
import scipy
from scipy import integrate
import math as m
import random as rd
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.pylab as pl

# Parameters for population of binaries circularized
nb = 10**5				# Choose 10^5 for my simulation, if you want to run a test, change it to ~10^2-10^3, *much* faster 
eta_val=2.
age_min = 1.				# Age bounds for my simulation, gives range of binary ages.  Want one value for your simulation
age_max = 10.

# Creating population of binaries to circularize
mu_pop = np.random.power(2.,size=nb)
A=1.75
B=2.01
e_pop = np.random.beta(A,B,size=nb)
Pmin=1.
Pmax=200.
xmin = np.log(Pmin)
xmax = np.log(Pmax)
lam = np.random.power(2.,size=nb)
x_pop = (1.-lam)*xmin + lam*xmax
P_pop = np.exp(x_pop)
age_pop = np.random.uniform(low=age_min, high=age_max, size=nb)

# The functions as defined in Leconte et al. (2010), A&A, doi:10.1051/0004-6361/201014337

def N(e):
	num = 1. + 15./2.*e**2. + 45./8.*e**4. + 5./16.*e**6.
	den = (1.-e**2.)**6.
	return num/den

def Na(e):
	num = 1. + 31./2.*e**2. + 255./8.*e**4. + 185./16.*e**6. + 25./64.*e**8.
	den = (1.-e**2.)**(15./2.)
	return num/den

def Omegae(e):
	num = 1. + 3./2.*e**2. + 1./8.*e**4.
	den = (1.-e**2.)**5.
	return num/den

def Ne(e):
	num = 1. + 15./4.*e**2. + 15./8.*e**4. + 5./64.*e**6.
	den = (1.-e**2.)**(13./2.)
	return num/den

def Omega(e):
	num = 1. + 3.*e**2. + 3./8.*e**4.
	den = (1.-e**2.)**(9./2.)
	return num/den

def Fe(e):
	om_eq = N(e)/Omega(e)			# This is the "pseudo-synchronous" rotation
	return Omegae(e)*om_eq - 18./11.*Ne(e)

def Fa(e):
	om_eq = N(e)/Omega(e)			# This is the "pseudo-synchronous" rotation
	return 4./11.*(N(e)*om_eq - Na(e))

# Ordinary differential equation integrate to get tidal evolution
# Parameters: 
# orb = (P[0],e[0]) initial orbital elements of binary
# time is time array you integrate simulation over
# mu is binary mass ratio M_2/M_1 *not* reduced mass
# eta is a powerlaw used for my problem, take out this parameter for your simulations

def TideEq(orb,time,mu,eta):
	P = orb[0]				# Orbital period, P = 2*pi/nb
	e = orb[1]				# Eccentricity
	tc = 0.3*(P/3.)**eta			# Circularization timescale, *need to change this for your problem*
	dedt = e*mu*(1.+mu)/tc*Fe(e)		# Differential equation for eccentricity
	dPdt = 1.5*P*mu*(1.+mu)/tc*Fa(e)	# Differential equation for orbital period
	return [dPdt,dedt]			# Returns differential equation for P and e

# *Note*, for your simulation, you likely want to change dPdt to dadt (semi-major axis evolution)

P_ev = []
e_ev = []
P_0 = []
e_0 = []
mu_ev = []
age_ev = []

for i in range(nb):
	# Individual parameters saving for circularized population of binaries
	Pb = P_pop[i]
	eb = e_pop[i]
	mub = mu_pop[i]
	ageb = age_pop[i]

	# This is where I tidally-circularize the population, you will likely want to do this for only one binary at a time	

	tb = np.linspace(age_min,ageb,num=3000)				# Array of times to calculate ODE over
	orbb = (Pb,eb)							# Initializing orbital elements
	solb = odeint(TideEq, orbb, tb, args=(mub,eta_val))		# Saving the time evolution of the arrays
	# Here, I save *only* the endpoint of the simulation.  You will likely want to keep the full time evolution
	if solb[-1,0]>0. and solb[-1,1]<1 and isinstance(solb[-1,0],float) and isinstance(solb[-1,1],float):
		P_ev.append(solb[-1,0])
		e_ev.append(solb[-1,1])
		P_0.append(solb[0,0])
		e_0.append(solb[0,1])
		mu_ev.append(mub)
		age_ev.append(ageb)

eta_str = r"eta={:.1f}".format(eta_val)
agemin_str = r"age0={:.1}".format(age_min)
data_save = (P_ev,e_ev,P_0,e_0,mu_ev,age_ev)
file_str = "Bin_Ev_ParTide_data_"+agemin_str+"_"+eta_str+".txt"
np.savetxt(file_str,data_save)

plt.hist(e_pop,density=True,bins=100,color='blue')
plt.hist(e_ev,density=True,bins=100,color='red')
plt.ylabel('Probability')
plt.xlabel('e')
	
plt.show()

plt.hist(np.log10(P_pop),density=True,bins=100,color='blue')
plt.hist(np.log10(P_ev),density=True,bins=100,color='red')
plt.ylabel('Probability')
plt.xlabel('P')
	
plt.show()

