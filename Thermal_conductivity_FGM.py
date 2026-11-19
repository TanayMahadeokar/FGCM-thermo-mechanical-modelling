import numpy as np
import matplotlib.pyplot as plt

# Parameters for the Matrix (Aluminium AA6082 - m)
E_m = 68.3 * 10**9          # Young's Modulus (Matrix) (Pa)
nu_m = 0.330                # Poisson's Ratio (Matrix) 
rho_m = 2.71             # Mass Density (Matrix) (g/cu. mm)
alpha_m = 23.1 * 10**-6  # Coefficient of linear thermal expansion (Matrix) (/K)
k_m = 172                # Thermal Conductivity (Matrix) (W/m-K)
Cv_m = 890               # Specific heat (Matrix) (J/kg-K)
K_m = 8.33 * 10**10       # Bulk Modulus (Matrix) (Pa)
G_m = 29 * 10**9         # Shear Modulus

# Parameters for the Reinforcement (PBT reinforced with 30% glass - c)
E_c = 10.48 * 10**9       # Young's Modulus (Reinforcement)
nu_c = 0.356            # Poisson's Ratio (Reinforcement)
rho_c = 1.55              # Mass Density (Reinforcement)
alpha_c = 32.5 * 10**-6  # Coefficient of linear thermal expansion (Reinforcement)
k_c = 0.245              # Thermal Conductivity (Reinforcement)
Cv_c = 1333             # Specific heat (Reinforcement)
K_c = 1.1 * 10**10       # Bulk Modulus (Reinforcement - Reinforced with 30% Glass Fibre)
G_c = 3.4955 * 10 ** 9    # Shear modulus

# Running variables used in the iteration
p = 0.40                 # Material parameter specifying that the FGM is metal-rich
w = 1                    # Width
X = [0]                  # Initialize X array
k = [k_c / (2 * k_m + k_c)]  # Initialize thermal conductivity array with initial value

# Defining the iteration equation for the thermal conductivity of the FGM
x_values = np.arange(0, w + 0.001, 0.001)  # Generate values from 0 to w with a step of 0.001
for x in x_values:
    V_m = (x / w)**p
    V_c = 1 - V_m
    # Calculating thermal conductivity for the current value of x
    k_x = k_c * ((1 + 3 * (k_m - k_c) * V_m) / (3 * k_c + (k_m - k_c) * V_c))
    k.append(k_x)
    X.append(x)

# Plotting the curve
plt.plot(X, k, label='Thermal Conductivity')
plt.xlabel('Distance from the PBT rich end (mm)')
plt.ylabel('Thermal Conductivity of the FGM (W/mK)')
plt.title('Thermal Conductivity of the Functionally Graded Material (FGM)')
plt.grid(True)
plt.legend()
plt.show()
