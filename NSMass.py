import numpy as np
import math
import matplotlib.pyplot as plt

# Constants

G = 6.6743E-11 # in N*m^2/Kg^2
c = 2.998E8 # in m/s
M_NS = (1.4)*1.989E30 # Approximate mass in Kg
R_NS = 11400 # Approximate radius in m

# Mass 
r = np.linspace(0, R_NS, num = 1001) # in m
    
M_BH = (r*(c**2))/(2*G) # in Kg
M_Tol = (M_NS)*(((5/2)*(r/(R_NS)*10**3)**3) - (3/2)*(r/(R_NS)*10**3)**5) # in Kg

outFile = open("PhantomBlackHoleMassData.txt", "w")
for i in range(0, 1001):
    outFile.write(str(M_BH[i]) + " " + str(M_Tol[i]) + " " + str(r[i]) + "\n")
outFile.close() 
