import numpy as np
import math

G = 6.6743E-11 # in N*m^2/Kg^2
c = 2.998E8 # in m/s
M_NS = 1.4 # in Solar Masses
R_NS = 10000 # Approximate radius in m
NDP = 1.44 # Pressure given in Solar Masses

r = np.linspace(0, 10000, num=1000) # in m
def P_Tol(r):
    return ((15*M_NS)/(8*math.pi*(R_NS)**3))*(1-(r/R_NS)**2)

#outFile = open("PhantomBlackHoleRadiusData.txt", "w")
#for i in range(0, 1000):
#    if (P_Tol[i] > NDP):
#        outFile.write("The required pressure exists within the Neutron Star at a radius of: " + 10*i + "\n")
#    else:
#        outFile.write("The required pressure does not exist within the Neutron Star." + "\n")
#outFile.close()

M = np.linspace(0, 3, num=1000) # in Solar Masses
def r_Sch(M):
    return (2*(M*(1.989E30))*G)/(c**2) # Includes conversion, outputs in m
def M_Tol(r_Sch):
    return M_NS*((5/2)*(r_Sch/R_NS)**3 - (3/2)*(r_Sch/R_NS)**5) # Outputs in Solar Masses

outFilenew = open("PhantomBlackHoleMassData.txt", "w")
for i in range(0, 1000):
    outFilenew.write(str(M[i]) + " " + str(M_Tol[i]) + "\n")
outFilenew.close() 