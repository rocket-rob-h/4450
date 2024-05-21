# convert from his mol table 31.59
#
# def shomate_cp(T, A, B, C, D, E):
#     return A + B * T + C * T**2 + D * T**3 + E / T**2
#
# # Shomate coefficients for OH
# A = 33.066178
# B = -11.363417
# C = 11.432816
# D = -2.772874
# E = -0.158558
#
# # Temperature in Kelvin
# T = 1200
#
# # Molar mass of OH in kg/mol
# MW_OH = 17 / 1000  # 17 g/mol
#
# # Calculate Cp for OH at 1200K in J/(mol·K)
# Cp_OH_mol = shomate_cp(T, A, B, C, D, E)
# print(f'Cp of OH at {T} K is {Cp_OH_mol:.3f} J/(mol·K)')
#
# # Convert Cp to J/(kg·K)
# OH_Cp = Cp_OH_mol / MW_OH
# print(f'Cp of OH at {T} K is {OH_Cp:.3f} J/(kg·K)')

#At 1200K
import numpy as np
H2_Cp = 15340
H2O_Cp = 2425
H_Cp = 20.79*1000 # per mole to kg as one g/mol
# OH_Cp = 2583.882
OH_Cp = 1858.235 # based on nist conversion



Ru = 8314 #J/kmol to leave molecular weight in g/mol. 1000 conversion done.

def calc_cp(Yi,Cpi):
    return Yi * Cpi
#
# Mole fractions from previous solve
# K_combined_3500: 0.1551675
# Progress variable xi: 0.1287019
# Mole Fractions:
XH2 = 0.5316159
XH2O = 0.4683841
XH = 0.0602819
XOH = 0.0602819
# Mass Fractions:
YH2 = 0.1011422
YH2O = 0.7963664
YH = 0.0057339
YOH = 0.0967575


MW_H2O = 18
MW_H = 1
MW_OH = 17
MW_H2 = 2.016

CpH2 = calc_cp(YH2,H2_Cp)
CpH = calc_cp(YH,H_Cp)
CpOH = calc_cp(YOH,OH_Cp)
CpH2O = calc_cp(YH2O,H2O_Cp)

print(f'CpH2 is {CpH2:.3f}')
print(f'CpH is {CpH:.3f}')
print(f'CpOH is {CpOH:.3f}')
print(f'CpH2O is {CpH2O:.3f}')

sumcpmix = CpH + CpOH + CpH2 + CpH2O

print(f'Cp of Mixture is {sumcpmix:.3f} J/kg.K')

summassfractions = YH + YOH + YH2 + YH2O

print(f'Mass fractions sum is {summassfractions:.3f} g/mol')


# Calculate the molar mass of the mixture
molar_mass_mixture = (XH2 * MW_H2 + XH2O * MW_H2O + XH * MW_H + XOH * MW_OH)
print(f'Molar Mass of Mixture is {molar_mass_mixture:.3f} g/mol')

def calc_R_specfic(Ru,molar_mass_mixture):
    return Ru / molar_mass_mixture

R_mixture = calc_R_specfic(Ru,molar_mass_mixture)
print(f'R_mixture is {R_mixture:.3f} J/kg.K')

def calc_gamma(Cp,R):
    return Cp / (Cp - R)

gamma_rocket = calc_gamma(sumcpmix,R_mixture)
print(f'The gamma for specific heats is {gamma_rocket:.4f}')







def sound_speed(gamma,R,Temperature):
    return np.sqrt(gamma * R * Temperature)

a_e = sound_speed(gamma_rocket,R_mixture,1200)
print(f'Exhaust sound speed is {a_e:.3f} m/s')

def Isentropic_relations_T(P1, T1, T2, gamma):
    """Calculate the ratio of temperatures in isentropic flow."""
    exp = (gamma/gamma-1)
    return P1 * ((T2/T1) ** exp)

Pe = 23802.760 # Exit Pressure/Inlet?
Po = 2e6 #  Stagnation pressure

# P2 = Isentropic_relations_T(Pe,1200,3500,gamma_rocket)
# print(f'P2 is {P2:.2f}Pa')


def calc_mach_number(P0, P, gamma):
    # Calculate the Mach number using the rearranged isentropic equation
    term = (P0 / P) ** ((gamma - 1) / gamma) - 1
    M = np.sqrt((2 / (gamma - 1)) * term)
    return M

M_e = calc_mach_number(Po,Pe,gamma_rocket)
print(f'Exit Mach number is {M_e:.3f}')

To = 3500 #stagnation temperature

def T1fromT2(To,gamma, Mach):
    return To / (1 + ((gamma - 1) / 2) * Mach**2)

T_e = T1fromT2(To,gamma_rocket,M_e)
print(f'T_e is {T_e:.3f} K')

def sound_speed(gamma,R,Temperature):
    return np.sqrt(gamma * R * Temperature)

a_e = sound_speed(gamma_rocket,R_mixture,T_e)
print(f'Exhaust sound speed is {a_e:.3f} m/s')

def Velocity(Mach, sound_speed):
    return Mach * sound_speed

V_e = Velocity(M_e,a_e)
print(f'Exit velocity is {V_e:.3f} m/s')

def Ideal_gas_law_rho(Pressure,R, Temperature):
    return Pressure / (R * Temperature)

rho_e = Ideal_gas_law_rho(Pe,R_mixture,T_e)
print(f'Density exit rho_e is {rho_e:.4e} kg/m^3')

