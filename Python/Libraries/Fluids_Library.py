#Preamble

import numpy as np
from scipy.optimize import fsolve, newton
# from scipy.optimize import newton
from math import *


#Initial Conditions
#Calculate Shock Pressures

# Driver Helium properties
gamma_4 = 1.667
R_4 = 2077  # J/kg.K
P_4 = 3e6  # 3Mpa
cp_4 = 5192.6  # J/kg.k
T_4 = 300  # K
L_4 = 0.23  # m
U_4 = 0

# Test Gas Air properties
gamma_1 = 1.4
R_1 = 287  # J/kg.K
P_1 = 15e3  # 15kpa
cp_1 = 1005
T_1 = 300  # K


def sound_speed(gamma,R,Temperature):
    return np.sqrt(gamma * R * Temperature)

def Ideal_gas_law_P(rho, R, Temperature):
    return rho * R * Temperature

def Ideal_gas_law_rho(Pressure,R, Temperature):
    return Pressure / (R * Temperature)

def Mach(velocity, sound_speed):
    return velocity / sound_speed

def Velocity(Mach, sound_speed):
    return Mach * sound_speed

#Shock Mach
def Shock_velocity(sound_speed, gamma, P2, P1):
    return sound_speed * np.sqrt ((gamma + 1) / (2 * gamma) * ((P2/P1)-1) + 1 )

def stagnation_pressure_ratio(gamma, Mach):
    return (1 + ((gamma - 1) / 2) * Mach ** 2) ** (gamma / (gamma - 1))

def stagnation_temperature_ratio(gamma, Mach):
    return (1 + ((gamma - 1) / 2) * Mach**2)

def Top_ratio(Denominator, Ratio):
    return Denominator * Ratio

def Bottom_ratio(Numerator, Ratio):
    return Numerator / Ratio

def Post_shock_velocity(sound_speed, gamma, P2, P1):
    numerator = (2 * gamma) / (gamma + 1)
    denominator = (P2 / P1) + ((gamma - 1) / (gamma + 1))
    return (sound_speed / gamma) * ((P2 / P1) - 1) * np.sqrt(numerator / denominator)

def Mach_reflected(M_2, gamma_1):
    first_term = M_2 * ((gamma_1 + 1)/4)
    second_term = np.sqrt((M_2**2 * (gamma_1 + 1) **2 ) / (16) + 1)
    return first_term + second_term

# From comp flow table

def calculate_mach_number_downstream(M1, gamma):
    """
    Calculate the Mach number downstream of the shock (M2) given the upstream Mach number (M1).
    """
    M2 = np.sqrt((1 + (gamma - 1)/2 * M1**2) / (gamma * M1**2 - (gamma - 1)/2))
    return M2

def shock_P2_P1(M1, gamma):
    """
    Calculate the pressure ratio across the shock (p2/p1).
    """
    pressure_ratio = 1 + 2*gamma/(gamma + 1) * (M1**2 - 1)
    return pressure_ratio

def shock_density_ratio(M1, gamma):
    """
    Calculate the density ratio across the shock (ρ2/ρ1).
    """
    density_ratio = ((gamma + 1) * M1**2) / (2 + (gamma - 1) * M1**2 )
    return density_ratio

def T2_T1_ratio(M1, gamma):
    """
    Calculate the temperature ratio across the shock (T2/T1).
    """
    # Using previously defined functions to get p2/p1 and ρ2/ρ1
    p_ratio = shock_P2_P1(M1, gamma)
    rho_ratio = shock_density_ratio(M1, gamma)
    temperature_ratio = p_ratio / rho_ratio
    return temperature_ratio

def calculate_total_pressure_ratio(M1, gamma):
    """
    Calculate the total pressure ratio across the shock (p02/p01).
    """
    p_ratio = shock_P2_P1(M1, gamma)
    total_pressure_ratio = (p_ratio * ((1 / (gamma * M1**2) + (gamma + 1) / 2) / ((gamma + 1) / (2 * gamma * M1**2))) ** (-gamma / (gamma - 1)))
    return total_pressure_ratio

def calculate_p1_p02(M1, gamma):
    """
    Calculate the total pressure ratio from state 1 to state 2 (p1/p02).
    """
    p02_p01 = shock_P2_P1(M1, gamma)
    p1_p01 = 1 / (1 + (gamma - 1) / 2 * M1**2) ** (gamma / (gamma - 1))
    p1_p02 = p1_p01 / p02_p01
    return p1_p02



# # Function to solve for P2/P1 using the equation provided
# def solve_for_P2_P1(x):
#     P2_P1 = x[0]
#     term_inside_brackets = (gamma_4 - 1) * (a_1 / a_4) * (P2_P1 - 1) / np.sqrt(2 * gamma_1 * (2 * gamma_1 + (gamma_1 + 1) * (P2_P1 - 1)))
#     return [P4_P1_known - P2_P1 * (1 - term_inside_brackets) ** (-2 * gamma_4 / (gamma_4 - 1))]
#
# # Initial guess for P2_P1
# P2_P1_initial_guess = [20.0]  # This is a guess and may need to be adjusted
# # Use fsolve to find the root
# P2_P1_solution = fsolve(solve_for_P2_P1, P2_P1_initial_guess)
#
# print(f"The calculated P2/P1 ratio is: {P2_P1_solution[0]:.4f}")
# # Rearrange the ratio to find Pressure of State 2 of accelerated test gas
# P_2 = P_1 * P2_P1_solution[0]



# OUTPUTS:
M_s =  4.14696203       # Upstream Mach number
Mb_2 =  0.43120352       # Downstream Mach number
Po2_Po1 =  0.12263917 # Total pressure ratio across the shock
P1_Po2 =  0.04422986  # Pressure ratio of state 1 to total pressure behind shock
rho2_rho1 =  4.64848391  # Density ratio across the shock
T2_T1 =  4.28028654    # Temperature ratio across the shock


exit_r = (83.9 / 2)/1000
throat_r = (7 / 2 )/ 1000

A_Astar = (np.pi * (exit_r ** 2))/ ((np.pi * (throat_r ** 2)))


def mach_A_Astar(A_Astar, gamma, mach_guess):
    """
    Find the Mach number for a given area ratio A/A* using the fsolve method.
    :param A_Astar: Known area ratio A/A*
    :param gamma: Specific heat ratio (γ)
    :param mach_guess: Initial guess for the Mach number
    :return: Mach number corresponding to the given area ratio
    """

    # Inner function to represent the equation to be solved
    def equation(M, gamma, A_Astar):
        term = (2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M ** 2)
        return (1 / M) * (term ** ((gamma + 1) / (2 * (gamma - 1)))) - A_Astar

    # Solve the equation
    mach_number, info, ier, msg = fsolve(equation, mach_guess, args=(gamma, A_Astar), full_output=True)
    if ier == 1:  # fsolve successfully found a root
        return mach_number[0]
    else:
        raise ValueError(f"Solution did not converge: {msg}")


# Find the Mach number for the given area ratio
# mach = mach_A_Astar(A_Astar_example, gamma, mach_guess)

###########################################################################3

def area_ratio(M, gamma, A_Astar):
    """ Function to calculate the deviation of calculated area ratio from the target. """
    return (1 / M) * ((2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M ** 2)) ** (
                (gamma + 1) / (2 * (gamma - 1)))) - A_Astar


def find_mach_number(A_Astar, gamma, mach_guess):
    """
    Find the Mach number for a given area ratio A/A* using the fsolve method.

    :param A_Astar: Known area ratio A/A*
    :param gamma: Specific heat ratio (γ)
    :param mach_guess: Initial guess for the Mach number
    :return: Mach number corresponding to the given area ratio
    """
    mach_number, info, ier, msg = fsolve(area_ratio, mach_guess, args=(gamma, A_Astar), full_output=True)
    if ier == 1:  # fsolve successfully found a root
        return mach_number[0]
    else:
        raise ValueError("Solution did not converge: " + msg)

gamma = 1.4  # Specific heat ratio for air
A_Astar = A_Astar # Example target area ratio A/A*
mach_guess = 5  # Starting guess for the Mach number

mach = find_mach_number(A_Astar, gamma, mach_guess)
# print(f"The Mach number for A/A* = {A_Astar} with γ = {gamma} is approximately {mach:.4f}")
#

def T5_T2wopress(gamma_2, M_2):  # Without pressure
    """
    Calculate the temperature ratio T5/T2 across a normal shock.

    Parameters:
    gamma_2: Specific heat ratio after the shock (γ).
    M_2: Mach number after the shock.

    Returns:
    T5/T2: Temperature ratio across the shock without considering pressure terms explicitly.
    """
    term_1 = (1 + ((2 * gamma_2) / (gamma_2 + 1)) * (M_2**2 - 1))
    term_2 = (2 + (gamma_2 - 1) * M_2**2) / ((gamma_2 + 1) * M_2**2)
    T5_T2 = term_1 * term_2
    return T5_T2

def T2_T1wopress(gamma_1, M_1): #without pressure
    term_1 = ((1+(2 * gamma_1) / gamma_1 + 1) * (M_1**2 - 1))
    term_2 = ((2 + (gamma_1 - 1) *M_1**2) / (gamma_1 + 1) * M_1**2)
    return term_1 * term_2

def T2_T1wpress(P2_P1, gamma_1):
    term_1 = P2_P1
    term_2 = ((1+(gamma_1 + 1 / gamma_1 - 1) + P2_P1) )
    term_3 = ((1+ (gamma_1 + 1)/ gamma_1 - 1) * P2_P1)
    return term_1 * (term_2 / term_3)

def T5_T1(T_5, T_1):
    return T_5 / T_1

def find_mach_from_pressure_ratio(pressure_ratio, gamma, initial_guess):
    """
    Find the upstream Mach number (M1) for a given pressure ratio (p2/p1) in a normal shock.

    :param pressure_ratio: The pressure ratio across the shock (p2/p1).
    :param gamma: Specific heat ratio (γ), for air it is typically 1.4.
    :param initial_guess: Initial guess for the Mach number.
    :return: The upstream Mach number (M1).
    """

    def f(M1):
        # Implicit equation for the upstream Mach number based on pressure ratio
        return 1 + (2 * gamma / (gamma + 1)) * (M1 ** 2 - 1) - pressure_ratio

    # Use the Newton-Raphson method to solve for M1
    M1 = newton(f, initial_guess)
    return M1

def speed_sound3(U4, U3, gamma4, a_4):
    return ((U4 - U3) * (gamma4 - 1) / 2) + a_4



def EQN(U_2, a_3, gamma_4, P5_P2):
    gamma = gamma_4 + 1
    gamma_ = gamma_4 - 1
    numer = ((gamma/gamma_) - 1) * (P5_P2 - 1)
    denom = (1 + (gamma / gamma_)) * (1 + ( gamma / gamma_) * P5_P2)
    return U_2 - a_3 * numer / np.sqrt(denom)

def normal_shock_relations(M1, gamma):
    """Calculate downstream Mach number M2 in a normal shock."""
    M2 = np.sqrt((1 + 0.5 * (gamma - 1) * M1**2) / (gamma * M1**2 - 0.5 * (gamma - 1)))
    return M2

def prandtl_meyer_function(M, gamma):
    """Calculate Prandtl-Meyer function for given Mach number M and specific heat ratio gamma."""
    return np.sqrt((gamma + 1) / (gamma - 1)) * np.arctan(np.sqrt((gamma - 1) * (M**2 - 1) / (gamma + 1))) - np.arctan(np.sqrt(M**2 - 1))

def p0_p(M, g=1.4):
    """
    Total to static pressure ratio for an isentropic flow.

    M: Mach number
    g: ratio of specific heats
    Returns: p0/p
    """
    return (T0_T(M, g))**( g / (g - 1.0) )

def calculate_EQN(u2, a3, gamma4, p5, p2):
    """
    Calculates the EQN based on the provided equation

    Parameters:
    u2 : float
        Velocity u2 from the equation.
    a3 : float
        Speed of sound a3 from the equation.
    gamma4 : float
        Ratio of specific heats gamma4 from the equation.
    p5 : float
        Pressure p5 from the equation.
    p2 : float
        Pressure p2 from the equation.

    Returns:
    float
        The calculated value of EQN.
    """
    numerator = ((gamma4 + 1) / (gamma4 - 1) - 1) * ((p5 / p2) - 1)
    denominator = np.sqrt((1 + (gamma4 + 1) / (gamma4 - 1)) * (1 + (gamma4 + 1) / (gamma4 - 1) * p5 / p2))
    EQN = u2 - a3 * (numerator / denominator)
    return EQN


# Oblique Shock Theory
# =============================================================================================================
# needs [0] list at end to print
def M1_n(M1, beta):
    return M1 * np.sin(beta)

def M2(M2_n, beta, theta):
    num = M2_n
    denom = np.sin(beta - theta)
    return num /denom


def beta_obl(M1, theta, g=1.4, tol=1.0e-6):
    sign_beta = np.sign(theta)
    theta = np.abs(theta)
    b1 = np.arcsin(1.0 / M1)
    b2 = b1 * 1.05

    def f_to_solve(beta, M1, theta, g):
        return theta_obl(M1, beta, g) - theta

    beta_solution = fsolve(f_to_solve, [b1, b2], args=(M1, theta, g), full_output=True, xtol=tol)[0]
    return sign_beta * beta_solution

def beta_obl2(M1, p2_p1, g=1.4):
    dum1 = np.sqrt(((g + 1.0) * p2_p1 + (g - 1.0)) / (2.0 * g))
    return np.arcsin(dum1 / M1)

def theta_obl(M1, beta, g=1.4):
    m1sb = M1 * np.sin(beta)
    t1 = (2.0 / np.tan(beta)) * (m1sb ** 2 - 1.0)
    t2 = (M1 ** 2) * (g + np.cos(2 * beta)) + 2.0
    return np.arctan2(t1, t2)

def M2_obl(M1, beta, theta, g=1.4):
    m1sb = M1 * np.sin(beta)
    numer = 1.0 + ((g - 1.0) / 2.0) * m1sb ** 2
    denom = g * m1sb ** 2 - (g - 1.0) / 2.0
    return np.sqrt(numer / (denom * (np.sin(beta - theta) ** 2)))

def rho2_rho1_obl(M1, beta, g=1.4):
    m1sb = M1 * np.sin(beta)
    numer = (g + 1.0) * m1sb ** 2
    denom = 2.0 + (g - 1.0) * m1sb ** 2
    return numer / denom

def p2_p1_obl(M1, beta, g=1.4):
    m1sb = M1 * np.sin(beta)
    return 1.0 + 2.0 * g / (g + 1.0) * (m1sb ** 2 - 1.0)

def T2_T1_obl(M1, beta, g=1.4):
    p2_p1_ratio = p2_p1_obl(M1, beta, g)
    rho2_rho1_ratio = rho2_rho1_obl(M1, beta, g)
    return p2_p1_ratio / rho2_rho1_ratio

def p02_p01_obl(M1, beta, g=1.4):
    m1sb = M1 * np.sin(beta)
    t1 = ((g + 1.0) / (2.0 * g * m1sb ** 2 - (g - 1.0))) ** (1.0 / (g - 1.0))
    t2 = (((g + 1.0) * m1sb ** 2) / (2.0 + (g - 1.0) * m1sb ** 2)) ** (g / (g - 1.0))
    return t1 * t2


############################
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

# # Example usage
# M1_example = 2.0
# theta_example = np.radians(30)
# beta_weak, beta_strong = solve_beta(M1_example, theta_example)
# print(f"Weak shock beta: {np.degrees(beta_weak):.2f} degrees")
# print(f"Strong shock beta: {np.degrees(beta_strong):.2f} degrees")