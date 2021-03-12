import numpy as np

G = 6.6743E-11 #in N*m^2/Kg^2
c = 2.998E8 #in m/s
M_NS = 1.4 #in Solar Masses
R_NS = 10000 #approximate radius in m

M = np.linspace(0, 3, 1000)

r_Sch = (2*(M*(1.989E30))*G)/(c**2)

M_Tol = M_NS*((5/2)*(r_Sch/R_NS)**3 - (3/2)*(r_Sch/R_NS)**5)

outFile = open("PhantomBlackHoleData.txt", "w")
outFile.write(str(M) + " " + str(M_Tol) + "\n")
outFile.close()
        