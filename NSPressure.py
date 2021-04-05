import numpy as np
import math
import matplotlib.pyplot as plt

c = 299792458 # speed of light in m/s
G = 6.6743E-11 # Newton's gravitational constant in m^3/kg*s^2
rho_central = 9.9E17 # central density in kg/m^3
R_NS = 11400 # Radius of neutron star in m
M_NS = 1.4*1.9891E30 # Mass of neutron star in kg

#Tolman VII Solution

r = np.linspace(0, R_NS, 1001)
xi = r/R_NS
rho_c = (15.0*M_NS)/(8.0*math.pi*R_NS**3)
C = ((8.0*math.pi/15.0)*R_NS**2*rho_c)*(G/c**2)

C_2 = math.atan(math.sqrt(C/(3.0*(1.0 - 2.0*C)))) + 0.5*math.log((1.0/6.0) + math.sqrt((1.0-2.0*C)/(3.0*C)))
E_lamb = 1.0 - C*xi**2*(5.0 - 3.0*xi**2)
phi_tol = C_2 - 0.5*np.log(xi**2 - (5.0/6.0) + np.sqrt(E_lamb/(3.0*C)))

p_tol = (1.0/(4.0*math.pi*R_NS**2))*((np.sqrt(3.0*C*E_lamb)*np.tan(phi_tol) - C*(5.0 - 3.0*xi**2)/2.0)/(4.0*math.pi*R_NS**2))

p_c = 0.9*p_tol[0]
print(p_c)

#Improved Tolman Solution

r = np.linspace(0, R_NS, 1001)

a_0 = 3.70625
a_1 = (-1.50266)
a_2 = 0.0643875
n = 0.903

alpha = a_0 + a_1*((C**n)/(rho_c*R_NS**2)) + a_2*((C**n)/(rho_c*R_NS**2))**2
C_2imp = math.atan(-(2.0*(10.0 - 3.0*alpha)*R_NS*math.sqrt(6.0*math.pi*rho_c*(15.0 - 16.0*math.pi*rho_c*R_NS**2)))/((48.0*math.pi*(10.0 - 3.0*alpha)*rho_c*R_NS) - 315.0)) + 0.5*math.log((1.0/6.0) + math.sqrt((5.0/(8.0*math.pi*rho_c*R_NS**2)) - (2.0/3.0)))
phi_imp = C_2imp - 0.5*np.log(xi**2 - (5.0/6.0) + np.sqrt((5.0*E_lamb)/(8*math.pi*rho_c*R_NS**2)))

p_imp = (np.sqrt(E_lamb_tol*rho_c/(10.0*math.pi))*(np.tan(phi_imp)/R_NS)) + ((1.0/15.0)*((3.0*xi**2) - 5.0)*rho_c) + ((6.0*(1.0-alpha)*rho_c)/((16.0*math.pi*(10.0 - (3.0*alpha))*rho_c*R_NS**2) - 105))

plt.plot(r/R_NS, p_tol/p_c)
plt.plot(r/R_NS, p_imp/p_c)
plt.show()