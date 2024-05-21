import numpy as np

# Interpolated values at 3500 K
ln_K_H2O_OH_H2 = -1.3320
ln_K_H2 = -1.0625

# Calculate ln(K_total)
ln_K_total = ln_K_H2O_OH_H2 + (ln_K_H2 / 2)

# Convert ln(K_total) to K_total
K_total = np.exp(ln_K_total)

# Constants
P_ref = 20  # bar = 2MPa
MW_H2O = 18.01528
MW_H = 1.00784
MW_OH = 17.007
MW_H2 = 2.01588
xi = 1.287019e-1  # Provided value

# Calculate mole fractions and mass fractions
def calculate_fractions(T, K, xi):
    X_H2 = 1.135 / 2.135
    X_H2O = 1 / 2.135
    X_H = xi / 2.135
    X_OH = xi / 2.135
    M = (X_H * MW_H) + (X_OH * MW_OH) + (X_H2O * MW_H2O) + (X_H2 * MW_H2)
    Y_H2 = (X_H2 * MW_H2) / M
    Y_H2O = (X_H2O * MW_H2O) / M
    Y_H = (X_H * MW_H) / M
    Y_OH = (X_OH * MW_OH) / M
    return {
        'Temperature': T,
        'K_combined': K,
        'xi': xi,
        'Mole Fractions': {
            'H2': X_H2,
            'H2O': X_H2O,
            'H': X_H,
            'OH': X_OH
        },
        'Mass Fractions': {
            'H2': Y_H2,
            'H2O': Y_H2O,
            'H': Y_H,
            'OH': Y_OH
        }
    }

# Calculate for T = 3500 K
T = 3500
result = calculate_fractions(T, K_total, xi)

# Print results
print(f"At {result['Temperature']} K, K_combined: {result['K_combined']:.7f}")
print(f"Progress variable xi: {result['xi']:.7f}")
print(f"Mole Fractions:")
for species, fraction in result['Mole Fractions'].items():
    print(f"  {species}: {fraction:.7f}")
print(f"Mass Fractions:")
for species, fraction in result['Mass Fractions'].items():
    print(f"  {species}: {fraction:.7f}")
