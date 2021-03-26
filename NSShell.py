import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import quadrature

G = 6.6743E-11 # in N*m^2/Kg^2
c = 2.998E8 # in m/s
M_NS = 2.3 # in Solar Masses
R_NS = 10000 # Approximate radius in m

M_BH = []
M_Tol = []
P_Grav = []
P_Degen = []

def Grav(r):
    return ((1/(r)**2)*(1-((r)/R_NS)**2))

# Given Radius,
# Find mass of a black hole of that Radius
# Find mass of Neutron Star with that Radius
# Compare

# Given radius,
# Find pressure due to gravity at that radius
# compare to Neutron Degeneracy Pressure.

r = np.linspace(0, 10000, num=1001) # in m
    
M_BH = (r*(c**2))/(2*G) # in Kg
M_Tol = (M_NS*1.989E30)*(((5/2)*(r/R_NS)**3) - (3/2)*(r/R_NS)**5) # in Kg

outFile = open("PhantomBlackHoleMassData.txt", "w")
for i in range(0, 1001):
    outFile.write(str(M_BH[i]) + " " + str(M_Tol[i]) + " " + str(r[i]) + "\n")
outFile.close() 

print("Check 1")

r = np.linspace(0, 10000, num=1001) # in m

#P_Grav, err = ((-(15*G*M_NS**2)/(8*math.pi*R_NS**3))*quadrature(Grav(r), 0, r)) # in N/m^2
P_Degen = ((((math.pi)**3)*(1.054E-34))/(15*(1.675E-27))*(((3*(M_NS*1.989E30))/(1.675E-27))/(((math.pi)**2)*(4/3)*(R_NS)**3))**(5/3)) # in N/m^2

outFilenew = open("PhantomBlackHolePressureData.txt", "w")
for i in range(0, 1001):
    if (i==0):
        P_Grav.append(0.0)
    else:
        result, err = quadrature(Grav, 0.001, r[i], tol=1e-6, rtol=1e-6)
        P_Grav.append(((-(15*G*M_NS**2)/(8*math.pi*R_NS**3))*result)) # in N/m^2
    outFilenew.write(str(P_Degen[i]) + " " + str(P_Grav[i]) + " " + str(r[i]) + "\n")
outFilenew.close()

Print("Check 2")