# Preamble
import numpy as np
import Fluids_Library
# import Propulsion_Library
from scipy.optimize import fsolve, newton

#1. Inlet Computations
'''a) Compute the angles at which the 3 shocks in the inlet occur relative to the flow'''

# Initial Conditions
Mach_inf = 15
T_inf = 220 # K
dynamic_pressure_q = 50000 # Pa
H2_Tot = 120e6 #J/kg Total Enthalpy 2
gamma_1 = 1.4 # Cp/Cv
R_air = 287 # J/kg.K

'''Assumptions
Adiabatic Walls, no heat transfer from flow to walls or vice-versa
No friction at wall inlets, two ramp intake
'''

theta_1_up = np.radians(5)
theta_2_up = np.radians(9) # 5 + 4 = 9
# 3rd shock reflects from cowl and turns flow back to parallel state 3
theta_3_down = np.radians(9)

def sound_speed(gamma,R,Temperature):
    return np.sqrt(gamma * R * Temperature)

# Sound speed to get a free stream velocity from the free stream temperature
print()

sound_speed_inf = sound_speed(gamma_1,R_air,T_inf)
print(f'a Free stream Sound Speed is {sound_speed_inf:.3f} m/s')

def Velocity(Mach, sound_speed):
    return Mach * sound_speed

velocity_inf = Velocity(Mach_inf,sound_speed_inf)
print(f'v Free stream Velocity is {velocity_inf:.3f} m/s')

def calc_density(dynamic_pressure, velocity):
    """Calculate the density from dynamic pressure and velocity."""
    return 2 * dynamic_pressure / (velocity ** 2)

rho_inf = calc_density(dynamic_pressure_q,velocity_inf)
print(f'rho The free stream density for this velocity and dynamic pressure is {rho_inf:.6f} kg/m^3')

def Ideal_gas_law_P(rho, R, Temperature):
    return rho * R * Temperature

P_static_inf = Ideal_gas_law_P(rho_inf,R_air,T_inf)
print(f'p Freestream Static Pressure from Ideal gas law is {P_static_inf:.3f} Pa')
def calc_dynamic_pressure(gamma, Pressure, Mach):
    return (gamma * Pressure * (Mach **2)) / 2

check_q = calc_dynamic_pressure(gamma_1, P_static_inf, Mach_inf)
print(f'Dynamic pressure is {check_q:.2f} Pa')
print(f'Therefore even with a higher static pressure than the US Standard atmosphere suggests for 220K, 23km atl, 0.005 kg, this is dynamic pressure driven')
print()
# Get all the values post oblique shock of 9 degrees at these initial conditions

def solve_beta(M1, theta, gamma=1.4):
    # Function to iterate
    def f_to_solve(beta, M1, theta, gamma):
        left_hand_side = np.tan(beta)
        right_hand_side = (((gamma + 1) * M1**2 * np.sin(beta)**2) /
                           (2 + (gamma - 1) * M1**2 * np.sin(beta)**2)) * np.tan(beta - theta)
        return left_hand_side - right_hand_side

    # Weak shock initial guess (close to theta)
    beta_weak_guess = theta + np.radians(1)  # A small angle above theta
    # Strong shock initial guess (further from theta)
    beta_strong_guess = np.pi/2 - np.radians(1)  # Near 90 degrees but within range

    # Solve for weak shock
    beta_weak_solution = fsolve(f_to_solve, beta_weak_guess, args=(M1, theta, gamma))[0]
    # Solve for strong shock
    beta_strong_solution = fsolve(f_to_solve, beta_strong_guess, args=(M1, theta, gamma))[0]

    return beta_weak_solution, beta_strong_solution

beta_weak, beta_strong = solve_beta(Mach_inf, theta_1_up)
print('--First Inlet Turn --'*3)
print()
print(f"Weak shock beta: {np.degrees(beta_weak):.2f} degrees")
print(f"Strong shock beta: {np.degrees(beta_strong):.2f} degrees")


def M2_obl(M1, beta, theta, g=1.4):
    m1sb = M1 * np.sin(beta)
    numer = 1.0 + ((g - 1.0) / 2.0) * m1sb ** 2
    denom = g * m1sb ** 2 - (g - 1.0) / 2.0
    return np.sqrt(numer / (denom * (np.sin(beta - theta) ** 2)))

M2_oblqiue_theta_1 = M2_obl(Mach_inf,beta_weak,theta_1_up,gamma_1)
print(f'The M2_oblqiue_theta_1 is {M2_oblqiue_theta_1:.2f}')

def M1_n(M1, beta):
    return M1 * np.sin(beta)

M1_n_inlet_1 = M1_n(Mach_inf,beta_weak)
print(f'The M1_n_inlet_1 Normal Mach 1 is {M1_n_inlet_1:2f}')


def ObliqueShock(M1, rho1, p1, theta, gamma_1=1.4):
    """
    Calculates properties behind an oblique shock wave given the upstream Mach number (M1),
    upstream density (rho1), upstream pressure (p1), deflection angle (theta in degrees),
    and specific heat ratio (gamma_1).
    """
    # Solve for weak beta using the solve_beta function
    beta_weak, _ = solve_beta(M1, theta, gamma_1)

    # Convert theta from degrees to radians
    # theta = np.radians(theta_deg)

    beta = beta_weak
    M1n = M1 * np.sin(beta)
    M1n2 = M1n * M1n
    gamma_1m1 = gamma_1 - 1.0
    gamma_1p1 = gamma_1 + 1.0

    M2n = np.sqrt((M1n2 + 2.0 / gamma_1m1) / (2.0 * gamma_1 / gamma_1m1 * M1n2 - 1))
    M2 = M2n / np.sin(beta - theta)
    M2t = M2 * np.cos(beta - theta)
    pratio = (2 * gamma_1 * M1n2 - gamma_1m1) / gamma_1p1
    rratio = gamma_1p1 * M1n2 / (gamma_1m1 * M1n2 + 2.0)
    p2 = p1 * pratio
    rho2 = rho1 * rratio
    # ptratio = pt1*sp.power(rratio,(k/km1))*sp.power(pratio,(-1.0/km1))
    a2 = np.sqrt(gamma_1 * p2 / rho2)

    return M2, M2n, M2t, rho2, p2, a2, beta

ObliqueShock_Inlet_1 = ObliqueShock(Mach_inf, rho_inf, P_static_inf, theta_1_up, gamma_1)

# Assign the results to individual variables
M2_oblique_theta1, M2n_oblique_theta1, M2t_oblique_theta1, rho2_oblique_theta1, P2_oblique_theta1, a2_oblique_theta1, beta_oblique_theta1 = ObliqueShock_Inlet_1

# Print the results with variable names
variables = ["M2", "M2n", "M2t", "rho2 (kg/m^3)", "p2 (Pa)", "a2 (m/s)", "beta (deg)"]
values = [M2_oblique_theta1, M2n_oblique_theta1, M2t_oblique_theta1, rho2_oblique_theta1, P2_oblique_theta1, a2_oblique_theta1, np.degrees(beta_oblique_theta1)]

for var_name, value in zip(variables, values):
    print(f"{var_name}: {value:.9f}")

# Function to calculate the temperature ratio across the oblique shock
def oblique_temp_ratio(p2, p1, rho2, rho1):
    temp_ratio = (p2 / p1) * (rho1 / rho2)
    return temp_ratio

# Calculate the temperature ratio across the first oblique shock
T2_T1_oblique_theta_1 = oblique_temp_ratio(P2_oblique_theta1, P_static_inf, rho2_oblique_theta1, rho_inf)
print(f'The temperature ratio across the first oblique shock is {T2_T1_oblique_theta_1:.2f}')

# The next variables for Inlet Turn 2
print('The next variables for Inlet Turn 2')

'''Now start the next set using the P2, T2 and M2 values for the next turn of 4 degrees'''