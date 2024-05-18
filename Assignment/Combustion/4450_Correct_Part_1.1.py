#Task 1.py

import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import numpy as np


T_o = 298 #K
Hr_O2 = 0
Hr_H2 = 0
Hfp_H2o_neg = -242000 # J/mol
Hfp_H2o_pos = 242000 #J/mol
Cp_h2o = 41.3 #J/(mol·K)
Cp_h2 = 30.2  #J/(mol·K)
T_limit = 3500 # K


def Ta(T_o, H_r, Hfp_H2o_neg, Cp_h2o):
    return T_o + ((H_r - Hfp_H2o_neg)/Cp_h2o)

# Adiabatic Stoichiometric ratio

Stoichyboi = Ta(T_o, (Hr_H2 + Hr_O2), Hfp_H2o_neg, Cp_h2o)

print(f"The adiabatic stoichiometric flame temperature for Hydrogen and Oxygen Combustion is {Stoichyboi:.2f} k")

# Now find the F/O ratio for when the flame temp drops below 3500 K

def find_phi(T_o, T_limit, Hfp_H2o_pos, Cp_h2o, Cp_h2):
    from scipy.optimize import fsolve

    def eq_solve(phi):
        x = 2 * phi - 2  # Excess hydrogen based on phi
        return T_o + (2 * Hfp_H2o_pos) / (2 * Cp_h2o + x * Cp_h2) - T_limit

    # Initial guess for phi (slightly above stoichiometric)
    initial_guess = 1.5
    phi_solution = fsolve(eq_solve, initial_guess)

    return phi_solution[0]

# Calculating the equivalence ratio
h2_rich_phi = find_phi(T_o, T_limit, Hfp_H2o_pos, Cp_h2o, Cp_h2)
print(f"Equivalence ratio at which the flame temperature reaches 3500 K is: {h2_rich_phi:.3f}")

def calculate_flame_temperature(T_o, Hfp_H2o_pos, Cp_h2o, Cp_h2, phi):
    x = 2 * phi - 2  # Excess hydrogen calculated from phi
    flame_temp = T_o + (2 * Hfp_H2o_pos) / (2 * Cp_h2o + x * Cp_h2)
    return flame_temp

# phi values from 1 to 10 as an example
phis = np.linspace(1, 10, 100)
flame_temps = [calculate_flame_temperature(T_o, Hfp_H2o_pos, Cp_h2o, Cp_h2, phi) for phi in phis]

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(phis, flame_temps, label='Adiabatic Flame Temperature')
plt.axvline(h2_rich_phi, color='red', linestyle='--', label='Eq ratio 2.135 at 3500K')
plt.scatter(1,Stoichyboi, color = 'blue', s = 50, label = '6157.56K at Eq = 1')
plt.text(1, Stoichyboi + 70, f'{Stoichyboi:.2f} K at $\phi$ = 1', color='blue')
plt.text(h2_rich_phi + 0.1, 3500 + 50, f'3500 K at $\lambda$ = {h2_rich_phi:.3f}', color='red')
plt.xlabel('Equivalence Ratio (λ)')
plt.ylabel('Flame Temperature (K)')
plt.legend()
plt.grid(True, axis ='y')
plt.show()


