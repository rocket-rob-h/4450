import numpy as np
import matplotlib.pyplot as plt
from K_Table import equilibrium_constants  # Library of K values

# Define the specific temperatures you're interested in
specific_temperatures = [1400, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3600]

K_combined = []

print("-- The Equilibrium constants have been found at each temperature -- ")

# Loop over the task temperatures
for T in specific_temperatures:
    if T in equilibrium_constants["H2O ⇋ OH + 1/2 H2"] and T in equilibrium_constants["H2 ⇋ 2H"]:
        ln_K_total = equilibrium_constants["H2O ⇋ OH + 1/2 H2"][T] + (equilibrium_constants["H2 ⇋ 2H"][T] / 2)
        K = np.exp(ln_K_total)
        K_combined.append((T, K))
        print(f"At {T} K, the combined K value is: {K:.4e}")

print()

# Constants
P_ref = 20  #bar = 2MPa
# MW_H2O = 18.01528
# MW_H = 1.00784
# MW_OH = 17.007
# M_H2 = 2.01588
MW_H2O = 18
MW_H = 1
MW_OH = 17
M_H2 = 2

# calculate mole fractions and mass fractions
def calculate_fractions(T, K):
    xi = np.sqrt(2 * K * (1 / P_ref))
    X_H2 = 0.5
    X_H2O = 0.5
    X_H = xi / 2
    X_OH = xi / 2
    M = (X_H * MW_H) + (X_OH * MW_OH) + (X_H2O * MW_H2O) + (X_H2 * M_H2)
    Y_H = (X_H * MW_H) / M
    Y_OH = (X_OH * MW_OH) / M
    Y_H_OH = Y_H + Y_OH
    return {
        'Temperature': T,
        'K_combined': K,
        'xi': xi,
        'Mole Fractions': {
            'H + OH': X_H + X_OH
        },
        'Mass Fractions': {
            'H + OH': Y_H_OH
        }
    }

# Calculate
results = [calculate_fractions(T, K) for T, K in K_combined]

# Print results
for result in results:
    print(f"At {result['Temperature']} K, K_combined: {result['K_combined']:.4e}")
    print(f"  xi: {result['xi']:.4e}")
    print(f"  Mole Fractions:")
    for species, fraction in result['Mole Fractions'].items():
        print(f"    {species}: {fraction:.4e}")
    print(f"  Mass Fractions:")
    for species, fraction in result['Mass Fractions'].items():
        print(f"    {species}: {fraction:.4e}")
    print()  #

# Plotting the results
temperatures, mass_fractions_H_OH = zip(*[(res['Temperature'], res['Mass Fractions']['H + OH']) for res in results])
plt.plot(temperatures, mass_fractions_H_OH, marker='o', label='Mass Fraction of H + OH')
plt.xlabel('Temperature (K)')
plt.ylabel('Mass Fraction of H + OH')
plt.title('Mass Fractions of H and OH vs. Temperature')
plt.grid(True)
plt.legend()
plt.savefig('MassFractions.png', dpi=600)
plt.show()
