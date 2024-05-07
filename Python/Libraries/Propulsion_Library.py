#This Library is used primarily for AERO 4450 at UQ towards the Final Assignment.
import numpy as np



#Equations from Final Exam Formula Sheet
def Ideal_P(rho, R, T):
    return rho * R * T

def gamma(Cp,Cv):
    return Cp/Cv

def Isentropic_relations_T(P1, T1, T2, gamma):
    exp = (gamma/gamma-1)
    return P1 * ((T2/T1) ** exp)

def Isentropic_relations_rho(P1,rho1, rho2, gamma):
    return P1 * (rho2/rho1)**gamma

def soundspeed(gamma, R, T):
    return np.sqrt(gamma * R * T)

def mach(v,a):
    return v * a

def stagnation_enthalpy(Cp, Tt):
    return Cp * Tt

def stagnation_temperature(gamma, M):
    return 1 + ((gamma-1)/2) * M**2

def stagnation_pressure(gamma, M):
    exp = (gamma / gamma - 1)
    return (1 + ((gamma-1)/2)*M**2) ** exp

# Normal Shock Waves

#Ratio of static density across a normal shock wave
def shock_rho2_rho1(gamma, M):
    m = M**2
    num = (gamma+1)*m
    denom = 2+(gamma-1)*m
    return num/denom

#Ratio of static pressure across a normal shock wave
def shock_P2_P1(gamma, M):
    gam = gamma+1
    m=M**2
    num = (2 * gamma * m) - gam
    return num / gam

def mach_downstream(gamma, M1):
    M2 = np.sqrt((1 + (gamma - 1)/2 * M1**2) / (gamma * M1**2 - (gamma - 1)/2))
    return M2

#Compressible Frictionless Duct with Heat Addition (Rayleigh Flow)

def Tt2(Tt1, q, Cp):
    return Tt1 + (q/Cp)

#Ratio of current to thermally choked Tt
def Tt_Ttstar(gamma, M):
    m = M**2
    num = (1+gamma)*m * (2+(gamma - 1)*m)
    denom = (1+gamma * m)**2
    return num / denom

#Rocket Nozzles
def Thrust(mdot, ve, Pe, Pinf, Ae):
    return mdot * ve + (Pe - Pinf)*Ae

def ideal_exit_velocity(gamma, R, Tc, Pe, Pc):
    exp = ((gamma-1)/gamma)

