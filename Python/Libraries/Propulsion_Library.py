#This Library is used primarily for AERO 4450 at UQ towards the Final Assignment.
import sys
import math
import numpy as np
import scipy as sp
import ThetaBetaM as tbm
import scipy.optimize as opt
from scipy.integrate import odeint
import matplotlib.pyplot as plt



"""

filename: ObliqueShock.py


 [M2, M2n, M2t, rho2, p2, a2, beta] = ObliqueShock(M1, rho1, p1, theta, k=1.4)

 Description:
 Computes the full flow state downstream of an oblique shock.

 Inputs:

   theta = turning angle of flow across the shock
   M1    = upstream Mach number
   rho1  = upstream density
   p1    = upstream pressure

 As this begins to seach for a solution at the low shock angle of theta+0.01 radians
 this will find the weak shock angle.

 Notes:
 gamma = 1.4 is assumed.

 Written by Brad Wheatley (built on Matlab code by vince Wheatley).

 """




def ObliqueShock(M1, rho1, p1, theta, k=1.4):
    beta = tbm.calcBeta(theta, M1, k)
    M1n = M1 * sp.sin(beta)
    M1n2 = M1n * M1n
    km1 = k - 1.0
    kp1 = k + 1.0

    M2n = sp.sqrt((M1n2 + 2.0 / km1) / (2.0 * k / km1 * M1n2 - 1))
    M2 = M2n / sp.sin(beta - theta)
    M2t = M2 * sp.cos(beta - theta)
    pratio = (2 * k * M1n2 - km1) / kp1
    rratio = kp1 * M1n2 / (km1 * M1n2 + 2.0)
    p2 = p1 * pratio
    rho2 = rho1 * rratio
    #    ptratio = pt1*sp.power(rratio,(k/km1))*sp.power(pratio,(-1.0/km1))
    a2 = sp.sqrt(k * p2 / rho2)

    return M2, M2n, M2t, rho2, p2, a2, beta


"""

Filename: ThetaBetaM.py

Description:
    A function to find the angle of an oblique shock to the upstream flow (in radians) given the upstream Mach number (M1) and the flow turning angle across the shock (in radians) (theta).

Reference:
    A.  NACA Report 1135, 'Equations, Tables and Charts for Compressible Flow'.

Inputs:
    theta - deflection angle of flow across the shock (radians)
    M1    - upstream Mach number
    beta  - angle of shock to the upstream flow (radians)
    k     - specific heat ratio (default is 1.4)

Methods:
    Reference A - Oblique shock formulae

Functions:
    calcBeta:
        Calculates the shock angle for given theta, M1 and k.
        Checks that the deflection angle is realistic for M1 and k.
    SonicLimitBeta:
        Calculates the shock wave angle (in radians) for sonic flow behind an oblique shock.
    MaxThetaBeta:
        Calculates the shock wave angle (in radians) for maximum stream deflection behind an oblique shock wave.
    calcTheta:
        Calculates deflection angle for an oblique shock with upstream Mach number M1 and shock angle beta.
    ThetaBetaM:
        Calculates the beta-Theta-M relationship for an oblique shock


External Calls:
    scipy.sin - sine of an angle (radians)
    scipy.tan - tan of an angle (radians)
    scipy.optimize.brentq - Find  the root of a function using Brent's method.
    sys.exit - exit the interpreter when conditions are met for failure


Written by Brad Wheatley 26/2/2013

"""


def calcBeta(theta, M1, k=1.4):
    """Calculate Beta given M1 and Theta. Check that Theta is realistic for the upstream Mach number first."""

    test_beta = MaxThetaBeta(M1, k)
    max_theta = calcTheta(M1, test_beta, k)

    if theta > max_theta:
        print("Deflection angle is too high. Maximum deflection for M=", M1, "is ", max_theta, "radians")
        sys.exit(1)

    beta = opt.brentq(lambda beta: ThetaBetaM(theta, M1, beta, k), 0.01, 2.0, maxiter=100)

    return beta


def SonicLimitBeta(M1, k=1.4):
    """Calculate the shock wave angle (in radians) for sonic flow behind an oblique shock. Refer to Ref A equation 167"""

    sin2beta = (k + 1) * M1 * M1 - (3.0 - k)
    sin2beta += sp.sqrt((k + 1) * ((k + 1) * sp.power(M1, 4) - 2 * (3 - k) * M1 * M1 + (k + 9)))
    sin2beta *= 1 / (4 * k * M1 * M1)
    sinbeta = sp.sqrt(sin2beta)

    if sinbeta <= 1.0:
        beta = sp.arcsin(sinbeta)
    else:
        print("Error. Sin(Theta) outside limits. Sin(Theta) = ", sinbeta)

    return beta


def MaxThetaBeta(M1, k=1.4):
    """Calculate the shock wave angle (in radians) for maximum stream deflection behind an oblique shock wave. Refer to Ref A equation 168."""

    sin2beta = (k + 1) * M1 * M1 - 4.0
    sin2beta += sp.sqrt((k + 1) * ((k + 1) * sp.power(M1, 4) + 8 * (k - 1) * M1 * M1 + 16))
    sin2beta *= 1 / (4 * k * M1 * M1)
    sinbeta = sp.sqrt(sin2beta)

    if sinbeta <= 1.0:
        beta = sp.arcsin(sinbeta)
    else:
        print("Error. Sin(Theta) outside limits. Sin(Theta) = ", sinbeta)

    return beta


def calcTheta(M1, beta, k=1.4):
    """Calculates deflection angle for an oblique shock with upstream Mach number M1 and shock angle beta. Refer to Ref A equation 139."""

    if beta <= 0.0:
        print("Shock angle must be positive.")
        sys.exit(1)
    if beta > sp.pi / 2.0:
        print("Shock angle must be less than pi/2 radians.")
        sys.exit(1)
    sinbeta = sp.sin(beta)
    cotbeta = 1.0 / sp.tan(beta)

    K1 = M1 * M1 * sinbeta * sinbeta

    tantheta = 2 * cotbeta * (K1 - 1.0) / (2 + M1 * M1 * (k + 1 - 2 * sinbeta * sinbeta))

    theta = sp.arctan(tantheta)

    return theta


def ThetaBetaM(theta, M1, beta, k=1.4):
    """Calculates the beta-Theta-M relationship for an oblique shock"""

    sinbeta = sp.sin(beta)
    K1 = M1 * M1 * sinbeta * sinbeta

    ans = sp.tan(beta - theta) / sp.tan(beta) - (2 + (k - 1) * K1) / ((k + 1) * K1)

    return ans


#Combustor Program

# Combustor entry state
# Lecture
# Tt3 = 3000
# p3 = 50.0e3
# M3 = 5
Tt3 = 3073.0
p3 = 28.74e3
M3 = 5.634
Cf = 0.002  # Skin friction factor in combustor

# Compute combustor area distribution
L = 0.5  # m, length of diverging combustor section

Nx = 50  # number of plot points
dx = L / Nx  # m, discretization size
x = np.linspace(0.0, L, Nx + 1)  # m, discretized combustor position
# A3 = 0.01              # m^2, initial combustor area
A3 = 0.0403  # m^2, initial combustor area

# Combustion chamber parameters
kb = 1.238
Rb = 290.3
cpb = kb * Rb / (kb - 1)

# Compute the remainder of the entry state
T3 = Tt3 / (1 + 0.5 * (kb - 1) * M3 ** 2)
Ht3 = cpb * Tt3
rho3 = p3 / (Rb * T3)
u3 = M3 * math.sqrt(kb * Rb * T3)
mdot = rho3 * u3 * A3

# 3-4: BURNER

print('3-4: COMBUSTOR')
km1o2 = 0.5 * (kb - 1.0)


# Define normalized stagnation temperature gradient
def dTtonTt(x, L, Tt3):
    X = x / L
    # Combustion chamber parameters
    hPR = 120.0e6  # for H2
    f = 0.0291  # stoiciometric H2
    phi = 0.75  # equivalence ratio
    etamax = 1.0  # max combustion effieciency
    theta = 1.0  # constant for combustion efficiency curve fit
    Tt = Tt3 + etamax * theta * X / (1 + (theta - 1) * X) * f * hPR * phi / cpb
    dTtdx = etamax * theta * (
                1 / (1 + (theta - 1) * X) - X * (theta - 1) / (1 + (theta - 1) * X) ** 2) * f * hPR * phi / (cpb * L)
    return Tt, dTtdx


# Define combustor area distribution
def A(x, A3, L):
    # return A3*(1.0+x/L)
    return A3 * (1.0 + 0.0 * x / L)


# Define normalized area gradient
def dAonA(x, A3, L):
    return A3 / L / A(x, A3, L)  # constant area combustor assumed


# Define effective diameter variation
def Deff(x, A3, L):
    return 2.0 * math.sqrt(A(x, A3, L) / math.pi)


# Define Mach number squared ODE (constant dAonA)
def dM2(M2, x, params):
    k, km1o2, Cf, dAonA, dTtonTt, Tt3, Deff, A3, L = params
    Tt, dTtdx = dTtonTt(x, L, Tt3)
    dM2 = M2 * (1.0 + km1o2 * M2) / (1.0 - M2) * (- 2.0 * dAonA(x, A3, L)
                                                  + (1.0 + k * M2) * dTtdx / Tt
                                                  + k * M2 * 4.0 * Cf / Deff(x, A3, L))
    return dM2


params = kb, km1o2, Cf, dAonA, dTtonTt, Tt3, Deff, A3, L
M2 = odeint(dM2, M3 * M3, x, args=(params,))
M = np.sqrt(M2)

# Compute other quantities
Tt, dTtdx = dTtonTt(x, L, Tt3)
Area = A(x, A3, L)
T = 0.0 * x.copy()
v = T.copy()
rho = T.copy()
p = T.copy()
for n in range(x.size):
    T[n] = Tt[n] / (1 + km1o2 * M2[n])
    v[n] = M[n] * np.sqrt(kb * Rb * T[n])
    rho[n] = mdot / (v[n] * Area[n])
    p[n] = rho[n] * Rb * T[n]

# Specify the exit state 4:
M4 = M[-1].copy()
p4 = p[-1].copy()
T4 = T[-1].copy()
rho4 = rho[-1].copy()
Tt4 = Tt[-1].copy()
v4 = v[-1].copy()
pt4 = p4 * (Tt4 / T4) ** (kb / (kb - 1))

plt.figure()
plt.plot(x, Tt / Tt3, label='$T_t/T_{t3}$')
plt.plot(x, T / T[0], label='$T/T_3$')
plt.plot(x, p / p[0], label='$p/p_3$')
plt.plot(x, M, label='$M$')
plt.plot(x, Area / Area[0], label='$A/A_3$')
plt.xlabel('$x$ [m]')
plt.legend()
# plt.show()
plt.savefig('examplecombustor.png', dpi=300)

print('Combustor Exit Conditions:')
print('p = ', p4, ' Pa, T = ', T4, ' K, M = ', M4)
print('=============================================================')




#Equations from Final Exam Formula Sheet
def Ideal_P(rho, R, T):
    """Calculate the pressure using the ideal gas law."""
    return rho * R * T

def gamma(Cp,Cv):
    return Cp/Cv

def Isentropic_relations_T(P1, T1, T2, gamma):
    """Calculate the ratio of temperatures in isentropic flow."""
    exp = (gamma/gamma-1)
    return P1 * ((T2/T1) ** exp)

def Isentropic_relations_rho(P1,rho1, rho2, gamma):
    return P1 * (rho2/rho1)**gamma

def soundspeed(gamma, R, T):
    """Calculate the speed of sound for an ideal gas."""
    return np.sqrt(gamma * R * T)

def mach(v,a):
    """Calculate the Mach number given velocity and speed of sound."""
    return v * a

def stagnation_enthalpy(Cp, Tt):
    """Calculate the stagnation enthalpy."""
    return Cp * Tt

def stagnation_temperature(gamma, M):
    """Calculate the stagnation temperature ratio."""
    return 1 + ((gamma-1)/2) * M**2

def stagnation_pressure(gamma, M):
    """Calculate the stagnation pressure ratio."""
    exp = (gamma / gamma - 1)
    return (1 + ((gamma-1)/2)*M**2) ** exp

# Normal Shock Waves

#Ratio of static density across a normal shock wave
def shock_rho2_rho1(gamma, M):
    """Calculate the ratio of static density across a normal shock."""
    m = M**2
    num = (gamma+1)*m
    denom = 2+(gamma-1)*m
    return num/denom

#Ratio of static pressure across a normal shock wave
def shock_P2_P1(gamma, M):
    """Calculate the ratio of static pressure across a normal shock."""
    gam = gamma+1
    m=M**2
    num = (2 * gamma * m) - gam
    return num / gam

def mach_downstream(gamma, M1):
    """Calculate the Mach number downstream of a normal shock wave."""
    M2 = np.sqrt((1 + (gamma - 1)/2 * M1**2) / (gamma * M1**2 - (gamma - 1)/2))
    return M2

def normal_shock_stagnation_pressure_ratio(M1, gamma):
    """Calculate the stagnation pressure ratio across a normal shock."""
    gam = gamma - 1
    m = M1**2
    num = (1 + (gamma - 1) / 2 * m)
    denom = (1 + gamma * (m - 1))
    return (num / denom)**(gamma / gam)

#Compressible Frictionless Duct with Heat Addition (Rayleigh Flow)

def Tt2(Tt1, q, Cp):
    return Tt1 + (q/Cp)

#Ratio of current to thermally choked Tt
def Tt_Ttstar(gamma, M):
    """Calculate the ratio of current to thermally choked stagnation temperatures in Rayleigh flow."""
    m = M**2
    num = (1+gamma)*m * (2+(gamma - 1)*m)
    denom = (1+gamma * m)**2
    return num / denom

#Rocket Nozzles
def Thrust(mdot, ve, Pe, Pinf, Ae):
    return mdot * ve + (Pe - Pinf)*Ae

def ideal_exit_velocity(gamma, R, Tc, Pe, Pc):
    exp = ((gamma-1)/gamma)

def Isp(F,mdot):
    return F / mdot * 9.81


############################################################################


def thrust_gas_turbines(m_dot_k, m_dot_f, v_e, v_0, p_e, p_inf, A_e):
    """Calculate thrust for gas turbines."""
    return m_dot_k * (v_e - v_0) + (p_e - p_inf) * A_e

def energy_balance(Q_dot, W_dot, m_dot_k, c_p, T_t_out, T_t_in):
    """Calculate the energy balance in a quasi-1D flow with heat addition."""
    return Q_dot - W_dot == m_dot_k * c_p * (T_t_out - T_t_in)

# def compressor_performance(pi_c, gamma, eta_c=None):
#     """Calculate non-ideal compressor performance measures."""
#     gam = gamma - 1
#     tau_c = pi_c ** (gam / gamma)
#     if eta_c is not None:
#         eta_c = (pi_c ** (gam / gamma) - 1) / (pi_c ** (eta_c * gam / gamma) - 1)
#     return tau_c, eta_c

def compressor_performance(pi_c, gamma, eta_c=None):
    """Calculate non-ideal compressor performance measures."""
    gam = gamma - 1
    tau_c = pi_c ** (gam / gamma)
    if eta_c is not None:
        eta_c = (pi_c ** (gam / gamma) - 1) / ((pi_c ** (eta_c * gam / gamma)) - 1)
    return tau_c, eta_c


# def turbine_performance(pi_t, gamma, eta_t=None):
#     """Calculate non-ideal turbine performance measures."""
#     gam = gamma - 1
#     tau_t = pi_t ** (-gam / gamma)
#     if eta_t is not None:
#         eta_t = 1 - (1 - pi_t ** ((1 - eta_t) * gam / gamma))
#     return tau_t, eta_t

def turbine_performance(pi_t, gamma, eta_t=None):
    """Calculate non-ideal turbine performance measures."""
    gam = gamma - 1
    tau_t = pi_t ** (-gam / gamma)
    if eta_t is not None:
        eta_t = 1 - (1 - pi_t ** ((1 - eta_t) * gam / gamma))
    return tau_t, eta_t

def burner_pressure_ratio(P4, P13, M3, M4, gamma):
    """Calculate burner pressure ratio."""
    numerator = 1 + gamma * M3**2 / 2
    denominator = 1 + gamma * M4**2 / 2
    return P4 / P13 * (numerator / denominator)

def burner_temperature_ratio_for_choking(M, gamma, tau_max):
    """Calculate the burner temperature ratio for choking in a ramjet."""
    return (1 + (gamma + 1) * M**2 / (2 + (gamma - 1) * M**2)) * tau_max

def maximum_equivalence_ratio(c_p, T_0, tau_b_max, f_st, H):
    """Calculate maximum equivalence ratio for a given burner temperature."""
    return c_p * T_0 * (tau_b_max - 1) / (f_st * H)

def a0(gamma, pi_b, pi_l, tau_b, tau_l, gamma, T_0, M_0):
    """Calculate speed of sound a_0 for specific thrust calculation."""
    term = (2 / (gamma - 1)) * ((pi_b * tau_b * pi_l * tau_l) ** (gamma / (gamma - 1)) * T_0 - 1)
    a_0_squared = (pi_b * tau_b * pi_l * tau_l) ** (-1 / gamma) * term - M_0
    return math.sqrt(a_0_squared)

def general_relation_burner_stagnation_temperature_ratio(c_p, T_t4, T_t3, f_st, H):
    """Calculate general relation between burner stagnation temperature ratio and equivalence ratio."""
    return c_p * (T_t4 - T_t3) / (f_st * H)


def burner_exit_mach_number(chi, gamma):
    """Calculate burner exit Mach number."""
    return 2 * chi / (1 - 2 * chi * gamma + (1 - 2 * chi * (gamma + 1)))

def specific_thrust_non_ideal(M_0, gamma, pi_b, pi_l, tau_b, tau_0):
    """Calculate specific thrust for non-ideal conditions."""
    term_inside_sqrt = (
            (2 / (gamma - 1)) *
            ((pi_b * tau_b) ** (gamma / (gamma - 1)) * tau_0 - 1)
    )

    a_0_squared_term = (pi_b * tau_b) ** (-1 / gamma) * term_inside_sqrt - M_0

    thrust = math.sqrt(a_0_squared_term)
    return thrust

############################################################################

def burner_stagnation_temperature_ratio(phi_s, H, cpb, T3b):
    """Calculate the burner stagnation temperature ratio."""
    tau_b = (phi_s * H) / (cpb * T3b) + 1
    return tau_b

def burner_exit_mach_number(chi, gamma):
    """Calculate the Mach number at the burner exit based on parameter chi and specific heat ratio gamma."""
    numerator = 2 * chi
    term_in_sqrt = 1 - 2 * chi * (gamma + 1)
    denominator = 1 - 2 * chi * gamma - math.sqrt(term_in_sqrt)
    M4 = math.sqrt(numerator / denominator)
    return M4

def burner_exit_chi(tau_b, M3b, gamma, Cf, Aw, A):
    """Calculate the parameter chi for the burner exit based on the given parameters."""
    numerator = tau_b * M3b**2 * (1 + (gamma - 1) / 2 * M3b**2)
    denominator_term = 1 - Cf * Aw / (24 * A)
    denominator = (1 + M3b**2 * denominator_term)**2
    chi = numerator / denominator
    return chi

def burner_exit_pressure(P3b, M3b, M4, gamma, Cf, Aw, A):
    """Calculate the burner exit pressure based on input parameters."""
    numerator_term = 1 - Cf * Aw / (24 * A)
    numerator = 1 + gamma * M3b**2 * numerator_term
    denominator = 1 + gamma * M4**2
    P4 = P3b * (numerator / denominator)
    return P4

def kinetic_energy_efficiency(M0, M3):
    """Calculate the kinetic energy efficiency based on the given Mach numbers at the inlet (M0) and at some point (M3)."""
    factor = (9 / M0)**0.7
    term = 0.018 * (1 - M3 / M0) + 0.12 * (1 - M3 / M0)**4
    eta_ke = 1 - factor * term
    return eta_ke

def inlet_pressure_recovery(psi, gamma, M0, M1):
    """Calculate inlet pressure recovery ratio."""
    pi_c = 1 - (1 - eta_c) * (1 + gamma * (M1 - M0)**2)
    return pi_c

def inlet_static_pressure_ratio(gamma, psi, M0, M3):
    """Calculate the static pressure ratio at the scramjet inlet."""
    pi_p = (1 + gamma * M0**2) / (1 + gamma * M3**2)
    return pi_p

def inlet_area_ratio(P0, M0, T3, P3, M3, T0, gamma, Cf, Aw, A):
    """Calculate the inlet area ratio."""
    T_ratio = T3 / T0
    pressure_term = P0 / P3
    mach_term = M0 / M3
    factor = (1 + gamma * M3**2 * (1 - Cf * Aw / (24 * A)))**2
    return pressure_term * mach_term * (T_ratio / factor)**0.5

def isentropic_nozzle_exit_mach_number(gamma, M4, P4, P10):
    """Calculate the isentropic nozzle exit Mach number using the given pressure and Mach number at one point."""
    # Calculate the pressure ratio part including the dynamic pressure component for M4
    pressure_ratio_component = (P4 / P10) * (1 + (gamma - 1) / 2 * M4 ** 2)
    M10_prime = math.sqrt((2 / (gamma - 1)) * (pressure_ratio_component - 1))
    return M10_prime

def isentropic_nozzle_exit_velocity(gamma, R, T4, M4):
    """Calculate the isentropic exit velocity of a nozzle."""
    V10_prime = (2 * gamma * R * T4 * (1 + ((gamma - 1) / 2) * M4**2))**0.5
    return V10_prime

def isentropic_nozzle_exit_temperature(T4, M4, M10_prime, gamma):
    """Calculate the isentropic nozzle exit temperature."""
    numerator = 1 + (gamma - 1) / 2 * M4**2
    denominator = 1 + (gamma - 1) / 2 * M10_prime**2
    T10_prime = T4 * (numerator / denominator)
    return T10_prime

###########################################################################################
#Assignment Equations

def growth_rate(ye, Ve, V3, ae, rho3, rhoe, r):
    """
    Calculate the growth rate of a mixing layer based on provided parameters.

    Parameters:
    ye (float): Specific heat ratio at exhaust.
    Ve (float): Velocity at the exhaust.
    V3 (float): Velocity at point 3.
    ae (float): Speed of sound at the exhaust.
    rho3 (float): Density at point 3.
    rhoe (float): Density at the exhaust.
    r (float): Density ratio parameter.

    Returns:
    float: Calculated growth rate delta prime.
    """
    Cd = 0.36  # Constant Cδ
    # Calculate π using the given formula
    pi = math.sqrt(ye - 1) * abs(Ve - V3) / ae

    # Calculate s using the given formula, inverting if necessary
    s = rho3 / rhoe
    if V3 > Ve:
        s = 1 / s  # Invert s if V3 > Ve as specified

    # Breakdown the formula into its components
    term1 = (1 + 4 * pi ** 2) ** -0.5  # (1 + 4π²)^-0.5
    term2 = (1 - r) / (1 + r) * (1 + math.sqrt(s))  # (1-r)/(1+r) * (1 + sqrt(s))
    term3 = (1 - math.sqrt(s)) / (1 + math.sqrt(s))  # (1-sqrt(s))/(1+sqrt(s))
    term4 = (1 + 2.9 * (1 + r)) / (1 - r)  # (1 + 2.9(1 + r))/(1 - r)

    # Calculate the growth rate delta prime
    delta_prime = term1 * Cd * s * term2 * term3 * term4

    return delta_prime


def assing_stagnation_enthalpy(mc, cpm, Tm, x, Lc, phi_r, mf, r, delta_prime, gamma_owl, gamma_cb, H):
    """
    Calculate the stagnation enthalpy of the total gas flow as a function of the distance downstream of the rocket exhaust.

    Parameters:
    mc (float): Mass flow rate of the combusted mixture.
    cpm (float): Specific heat at constant pressure of the mixture.
    Tm (float): Temperature of the mixture.
    x (float): Distance downstream from the rocket exhaust.
    Lc (float): Characteristic length.
    phi_r (float): Combustion efficiency.
    mf (float): Mass flow rate of the fuel.
    r (float): Fuel/air ratio.
    delta_prime (float): Growth rate parameter.
    gamma_owl (float): Specific heat ratio at the owl condition.
    gamma_cb (float): Specific heat ratio at the cb condition.
    H (float): Enthalpy per unit mass.

    Returns:
    float: Stagnation enthalpy at distance x.
    """
    # Term from the combusted mixture
    enthalpy_term = mc * cpm * Tm

    # Dynamic enthalpy term based on distance and other factors
    dynamic_term = (x / Lc + (phi_r - 1) * mf * r * delta_prime * x / (2 * (gamma_owl - gamma_cb))) * H

    # Total stagnation enthalpy
    stagnation_enthalpy = enthalpy_term + dynamic_term
    return stagnation_enthalpy

