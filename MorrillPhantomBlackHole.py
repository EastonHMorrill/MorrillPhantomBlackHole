import numpy as np
import math
import matplotlib.pyplot as plt

G = 6.6743E-11 # in N*m^2/Kg^2
c = 2.998E8 # in m/s
M_NS = 2.3 # in Solar Masses
R_NS = 10000 # Approximate radius in m

NDPCalc = (((math.pi)^3)*(1.054E-34))/(15*(1.675E-27))*(((3*(M_NS*1.989E30))/(1.675E-27))/(((math.pi)**2)*(4/3)*(R_NS)**3))**(5/3)

NDP = []
for i in range(0,1001):
    NDP.append(NDPCalc) # Pressure given in N/m^2

rho_Tol = []
P = []
r = np.linspace(0, 10000, num=1001) # in m
for i in range(0, 1001):
    if (r[i]==0):
        rho_Tol.append(0.0)
        P.append(0.0)
    else:
        rho_Tol.append(((15*M_NS)/(8*math.pi*(R_NS)**3))*(1-((r[i])/R_NS)**2)) # Output in Solar masses/m^3
        P.append((rho_Tol[i])*(G*M_NS/((r[i])**2))*((1.989E30)**2)*(R_NS-(r[i]))) # Output in N/m^2

outFile = open("PhantomBlackHolePressureData.txt", "w")
for i in range(0, 1001):
    if (P[i] > NDP[i]):
        outFile.write("The required pressure exists within the Neutron Star at a radius of: " + str(10*i) + "m" + "\n")
    else:
        outFile.write("The requred pressure does not exist at a radius of: " + str(10*i) + "m" + "\n")
outFile.close()

plt.plot(r, P)
plt.plot(r, NDP)
plt.savefig("PressurePlot.pdf")
plt.close()

M = np.linspace(0, 3, num=1001) # in Solar Masses
r_Sch = (2*(M*(1.989E30))*G)/(c**2) # Includes conversion, outputs in m
M_Tol = M_NS*((5/2)*(r_Sch/R_NS)**3 - (3/2)*(r_Sch/R_NS)**5) # Outputs in Solar Masses

outFilenew = open("PhantomBlackHoleMassData.txt", "w")
for i in range(0, 1001):
    outFilenew.write(str(M[i]) + " " + str(M_Tol[i]) + " " + str(r[i]) + "\n")
outFilenew.close() 

plt.plot(r_Sch, M)
plt.plot(r_Sch, M_Tol)
plt.savefig("MassPlot.pdf")
plt.close()