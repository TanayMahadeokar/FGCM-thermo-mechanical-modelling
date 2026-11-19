# THis plot contains the iteration engines of the following 4 models for thr FG-PBTRG
# Rule of Mixtures
# Modified Rule of Mixtures
# Mori-Tanaka 
# Haphin-Tsai

import matplotlib.pyplot as plt
import numpy as np

# Parameters for the Matrix (Aluminium AA6082 - m)
E_m = 68.3e9                # Young's Modulus (Matrix) (Pa)
nu_m = 0.330                # Poisson's Ratio (Matrix) 
rho_m = 2.71                # Mass Density (Matrix) (g/cu. mm)
alpha_m = 23.1e-6           # Coefficient of linear thermal expansion (Matrix) (/K)
k_m = 172                   # Thermal Conductivity (Matrix) (W/m-K)
Cv_m = 890                  # Specific heat (Matrix) (J/kg-K)
K_m = 8.33e10               # Bulk Modulus (Matrix) (Pa)
G_m = 29e9                  # Shear Modulus

# Parameters for the Reinforcement (PBT reinforced with 30% glass - c)
E_c = 10.48e9               # Young's Modulus (Reinforcement)
nu_c = 0.356                # Poisson's Ratio (Reinforcement)
rho_c = 1.55                # Mass Density (Reinforcement)
alpha_c = 32.5e-6           # Coefficient of linear thermal expansion (Reinforcement)
k_c = 0.245                 # Thermal Conductivity (Reinforcement)
Cv_c = 1333                 # Specific heat (Reinforcement)
K_c = 1.1e10                # Bulk Modulus (Reinforcement - Reinforced with 30% Glass Fibre)
G_c = 3.4955e9              # Shear modulus

# Running variables used in the iteration
p = 0.4                     # Material parameter specifying that the FGM is metal-rich
t = 1                       # thickness of the FGM (arbitrary units) (In our Model the thickness is 10 cm, scaled to 1)
l = 10                      # length of the FGM (arbitrary units) (In our Model the length is 100cm, scaled to 10) 
w = 5                       # width of the FGM (arbitrary units)  (In our model the width is 50 cm, scaled to 5)

# Defining the aspect ratio of fillers for the MROM
s = 10                      # s = length/thickness (arbitrary units)

# Defining the parameters for the Mori-Tanaka Model
K_m1 = E_m / (3 * (1 - 2 * nu_m))
K_c1 = E_c / (3 * (1 - 2 * nu_c))
G_m1 = E_m / (2 * (1 + nu_m))
G_c1 = E_c / (2 * (1 + nu_c))
f_m = G_m1 * (9 * K_m1 + 8 * G_m1) / (6 * (K_m1 + 2 * G_m1))

#Defining the paramters for the Halphin-Tsai model
zeta_11 = 0.60      # 2*(t/l) because we have modeled the Young's Modulus as a function of the width
zeta_22 = 0.25      # 2*(t/w) --------------------------------||------------------------------------
eta_11 = (((E_c / E_m) - 1) / (E_c / E_m) + zeta_11)
eta_22 = (((E_c / E_m) - 1) / (E_c / E_m) + zeta_22)
eta_12 = (((G_c / G_m) - 1) / (G_c / G_m))

# Defining the parameters for the Modified halphin-Tsai model

X1 = [0]                                                                                                                    # Initialising the X variable for ROM
X2 = [0]                                                                                                                    # Initialising the X variable for MROM
X3 = [0]                                                                                                                    # Initialising the X variable for the Mori-Tanaka Model
X4 = [0]                                                                                                                    # Initialising the X variable for Halpin-Tsai Model

E1 = [E_c]                                                                                                                  # Initial Young's Modulus for ROM
E2 = [25.3043 * E_c]                                                                                                        # Initial Young's Modulus for MROM
E3 = [(9 * K_m1 * G_m1) / (3 * K_m1 + G_m1)]                                                                                    # Initial Young's Modulus for the Mori-Tanaka Model
E4 = [(3 / 8) * ((1 + (zeta_11 * eta_11)) / (1 - eta_11)) * E_m + (5 / 8) * ((1 + zeta_22 * eta_22) / (1 - eta_22)) * E_m]  # Initial Young's Modulus for the Halpin-Tsai Model

# Defining the iteration equation for the Young's Modulus for the FGM using the Rule of Mixtures
x1_values = np.arange(0, t, 0.001) # Values from 0.001 to t with step size 0.001 
for x1 in x1_values:
    V_c1 = (x1 / t) ** p
    V_m1 = 1 - V_c1
    E_x1 = E_c * ((1 + ((E_m / E_c) - 1) * (V_m1) ** (2 / 3)) / (1 + ((E_m / E_c) - 1) * ((V_m1) ** (2 / 3) - V_m1)))
    E1.append(E_x1)
    X1.append(x1)

# Defining the iteration Equation for Young's Modulus using the Modified Rule of Mixtures (MROM) for the FGM
x2_values = np.arange(0, t, 0.001) # Values from 0.001 to w with step size 0.001
for x2 in x2_values:
    V_c2 = (x2 / t) ** p
    V_m2 = 1 - V_c2
    eta_zero = np.cos((np.pi/ (3 * w)) * x2) ** 4
    beta = np.sqrt((2 * E_m) / (E_c * (1 + nu_m) * np.log(1 / V_c2)))
    eta_one = 1 - ((np.tanh(beta * s)) / (beta * s)) 
    E_x2 = (eta_zero * eta_one * E_c * V_c2) + E_m * V_m2
    E2.append(E_x2)
    X2.append(x2)

# Defining the iteration Equation for Young's Modulus using the Mori-Tanaka Method for the FGM
x3_values = np.arange(0, t, 0.001) # Values from 0.001 to t with step size 0.001
for x3 in x3_values:
    V_c3 = (x3 / t) ** p
    V_m3 = 1 - V_c3 
    K_x3 = K_m1 + ((V_c3 * (K_c1 - K_m1)) / (1 + (1 - V_c3) * (3 * (K_c1 - K_m1) / (3 * K_m1 + 4 * G_m1))))
    G_x3 = G_m1 + ((V_c3 * (G_c1 - G_m1)) / (1 + (1 - V_c3) * ((G_c1 - G_m1) / (G_c1 + f_m))))
    E_x3 = ((9 * K_x3 * G_x3) / (3 * K_x3 + G_x3))
    E3.append(E_x3)
    X3.append(x3)

#Defining the iteration eqaution for Young's Modulus using the Halpin-Tsai Method for the FGM
x4_values = np.arange(0, t, 0.001) # Values from 0.001 to t with step size 0.001
for x4 in x4_values:
    V_c4 = (x4 / t) ** p
    V_m4 = 1 - V_c4
    E_transverse = 3 * (1 + zeta_11 * eta_11 * V_c4) / (1 - eta_11 * V_c4)
    E_longitudinal = 5 * (1 + zeta_22 * eta_22 * V_c4) / (1 - eta_22 * V_c4)
    E_x4 = (E_m / 8) * (E_transverse + E_longitudinal)
    E4.append(E_x4)
    X4.append(x4)
 
# Plotting the curves
plt.plot(X1, E1, label = 'Variation according to Rule of Mixtures')
plt.plot(X2, E2, label = 'Variation according to Modified Rule of Mixtures')
plt.plot(X3, E3, label = 'Variation according to Mori-Tanaka Model')
plt.plot(X4, E4, label = 'Variation according to Haphin-Tsai Model')
plt.xlim(-0.25, 1.25)
plt.ylim(bottom = -0.4e11, top = 0.8e11)
plt.xlabel('Distance from PBT rich end (mm)')
plt.ylabel("Young's Modulus of the FGM (Pa)")
plt.grid(True)
plt.legend()
plt.show() 