from propagator import Sf
from dirac import g5, Pp, gamma, Pm
from math import sqrt, log
import numpy as np
############################################################
#  Calculation of the bare mass
############################################################

def mb(x, z, L):
    """Calculate bare mass given with L,z,x as input."""
    return 1 - sqrt(- 2. * x* z / L + 1)

############################################################
#  CORRELATION FUNCTIONS
############################################################

def fX(L, T, x0, m, gam, thetah, thetal, gamp = g5):
    """Template for heavy light correlation functions a la fA, fP, kV
    etc."""
    pl = np.array([float(thetal)/L,]*3)
    ph = np.array([float(thetah)/L,]*3)
    Sh = lambda x, y : Sf(ph, x, y, T, m)
    Sl = lambda x, y : Sf(pl, x, y, T, 0)
    return 1.5 * np.trace( Sh(x0, 1) *
                           Pp * gamp * Sl(1, x0) * gam).real

def fA(L, T, x0, m, thetah, thetal):
    b_coeff = (3./2 - 0.5*sqrt(1 - 2.* m))
    return fX(L, T, x0, m, gamma[0] * g5, thetah, thetal)*b_coeff

def fP(L, T, x0, m, thetah, thetal):
    b_coeff = (3./2 - 0.5*sqrt(1 - 2.* m))
    return fX(L, T, x0, m, g5, thetah, thetal)*b_coeff

def kv(L, T, x0, m, thetah, thetal):
    b_coeff = (3./2 - 0.5*sqrt(1 - 2.* m))
    return fX(L, T, x0, m, gamma[1], thetah, thetal, gamma[1])*b_coeff

def f1(L, T, m, thetah, thetal):
    """f1, heavy light version"""
    pl = np.array([float(thetal)/L,]*3)
    ph = np.array([float(thetah)/L,]*3)
    Sh = lambda x, y : Sf(ph, x, y, T, m)
    Sl = lambda x, y : Sf(pl, x, y, T, 0)
    return 1.5 * np.trace( Sl(1, T-1) * g5 * Pp * 
                           Sh(T-1, 1) * g5 * Pm).real

def fXstat(L, T, x0, gam, thetal, gamp = g5):
    """Template for heavy light correlation functions a la fA, fP, kV
    etc."""
    pl = np.array([float(thetal)/L,]*3)
    Sh = lambda x, y : Pp
    Sl = lambda x, y : Sf(pl, x, y, T, 0)
    return 1.5 * np.trace( Sh(x0, 1) *
                           Pp * gamp * Sl(1, x0) * gam).real

def fAstat(L, T, x0, thetal):
    return fXstat(L, T, x0, gamma[0] * g5, thetal)

def f1stat(L, T, thetal):
    """f1, static light version"""
    pl = np.array([float(thetal)/L,]*3)
    Sh = lambda x, y : Pp
    Sl = lambda x, y : Sf(pl, x, y, T, 0)
    return 1.5 * np.trace( Sl(1, T-1) * g5 * Pp * 
                           Sh(T-1, 1) * g5 * Pm).real

def Okin(theta, L, dx):
    l = np.exp(1J * theta/L)
    return 3.*(l + 1./l - 2)*dx

def f1kin(L, T, thetah, thetal):
    pl = np.array([float(thetal)/L,]*3)
    Sh = lambda x, y : Okin(thetah, L, abs(x-y) + 1)*Pp
    Sl = lambda x, y : Sf(pl, x, y, T, 0)
    return 1.5 * np.trace( Sl(1, T-1) * g5 * Pp * 
                           Sh(T-1, 1) * g5 * Pm).real

def fXkin(L, T, x0, gam, thetah, thetal, gamp = g5):
    """Template for heavy light correlation functions a la fA, fP, kV
    etc."""
    pl = np.array([float(thetal)/L,]*3)
    Sh = lambda x, y : Pp*Okin(thetah, L, abs(x-y) + 1)
    Sl = lambda x, y : Sf(pl, x, y, T, 0)
    return 1.5 * np.trace( Sh(x0, 1) *
                           Pp * gamp * Sl(1, x0) * gam).real

def fAkin(L, T, x0, thetah, thetal):
    return fXkin(L, T, x0, gamma[0] * g5, thetah, thetal)

############################################################
#  Observables
############################################################

def R1(L, x, z, theta1, theta2):
    """Correlation function f1(theta1)/f1(theta2),
    T = L"""
    m = mb(x, z, L)
    ff1 = lambda theta : f1(L, L/2, m, theta, theta)
    return ff1(theta1) / ff1(theta2)

def logR1(L, x, z, theta1, theta2):
    """Correlation function R_1 = log(f1(theta1)/f1(theta2)),
    T = L"""
    m = mb(x, z, L)
    ff1 = lambda theta : f1(L, L/2, m, theta, theta)
    return log( ff1(theta1) / ff1(theta2) )

def RA(L, x, z, theta1, theta2):
    """Correlation function fA(theta1)/fA(theta2)"""
    m = mb(x, z, L)
    ffa = lambda theta : fA(L, L, L/2, m, theta, theta)
    return ffa(theta1) / ffa(theta2)

def logRA(L, x, z, theta1, theta2):
    """Correlation function R_A = log(fA(theta1)/fA(theta2))"""
    m = mb(x, z, L)
    ffa = lambda theta : fA(L, L, L/2, m, theta, theta)
    return log( ffa(theta1) / ffa(theta2) )

def RP(L, x, z, theta1, theta2):
    """Correlation function fP(theta1)/fP(theta2)"""
    m = mb(x, z, L)
    ffp = lambda theta : fP(L, L, L/2, m, theta, theta)
    return ffp(theta1) / ffp(theta2)

def logRP(L, x, z, theta1, theta2):
    """Correlation function R_P = log(fP(theta1)/fP(theta2))"""
    m = mb(x, z, L)
    ffp = lambda theta : fP(L, L, L/2, m, theta, theta)
    return log( ffp(theta1) / ffp(theta2) )

def RV(L, x, z, theta1, theta2):
    """Correlation function kV(theta1)/kV(theta2)"""
    m = mb(x, z, L)
    fkv = lambda theta : kv(L, L, L/2, m, theta, theta)
    return fkv(theta1) / fkv(theta2)

def logRV(L, x, z, theta1, theta2):
    """Correlation function R_V = log(kV(theta1)/kV(theta2))"""
    m = mb(x, z, L)
    fkv = lambda theta : kv(L, L, L/2, m, theta, theta)
    return log( fkv(theta1) / fkv(theta2) )

def LGammaPS(L, x, z, theta):
    """Correlation function L*Gamma_PS ~ (d/dx_0) log(fA(x_0))"""
    m = mb(x, z, L)
    return 0.5*L*(log(-fA(L, L, L/2 + 1, m, theta, theta)) \
           - log(-fA(L, L, L/2 - 1, m, theta, theta)))

def LGammaP(L, x, z, theta):
    """Correlation function L*Gamma_P ~ (d/dx_0) log(fP(x_0)), x_0=T/2"""
    m = mb(x, z, L)
    return 0.5*L*(log(fP(L, L, L/2 + 1, m, theta, theta)) \
           - log(fP(L, L, L/2 - 1, m, theta, theta)))

def LGammaV(L, x, z, theta):
    """Correlation function L*Gamma_V ~ (d/dx_0) log(kV(x_0)), x_0=T/2"""
    m = mb(x, z, L)
    return 0.5*L*(log(kv(L, L, L/2 + 1, m, theta, theta)) \
           - log(kv(L, L, L/2 - 1, m, theta, theta)))
    
def RPSP(L, x, z, theta):
    """Correlation function ~ fA(T/2)/fP(T/2)"""
    m = mb(x, z, L)
    return -fA(L, L, L/2, m, theta, theta)/ \
            fP(L, L, L/2, m, theta, theta)

def RPSV(L, x, z, theta):
    """Correlation function ~ fA(T/2)/kV(T/2)"""
    m = mb(x, z, L)
    return -fA(L, L, L/2, m, theta, theta)/ \
            kv(L, L, L/2, m, theta, theta)
    
def YPS(L, x, z, theta):
    """Correlation function Y_PS ~ fA(T/2)/sqrt(f1)"""
    m = mb(x, z, L)
    return -fA(L, L, L/2, m, theta, theta)/ \
            sqrt(f1(L, L, m, theta, theta))

def logYPS(L, x, z, theta):
    """Correlation function log(Y_PS) ~ log(fA(T/2)/sqrt(f1))"""
    m = mb(x, z, L)
    return log(-fA(L, L, L/2, m, theta, theta)/ \
                sqrt(f1(L, L, m, theta, theta)))

def YV(L, x, z, theta):
    """Correlation function Y_V ~ kV(T/2)/sqrt(k1)"""
    m = mb(x, z, L)
    return kv(L, L, L/2, m, theta, theta)/ \
           sqrt(f1(L, L, m, theta, theta))

def logYV(L, x, z, theta):
    """Correlation function log(Y_V) ~ log(kV(T/2)/sqrt(k1))"""
    m = mb(x, z, L)
    return log(kv(L, L, L/2, m, theta, theta)/ \
               sqrt(f1(L, L, m, theta, theta)))

def RYPSV(L, x, z, theta):
    """Correlation function Y_PS/Y_V"""
    m = mb(x, z, L)
    ffa = lambda theta : -fA(L, L, L/2, m, theta, theta)/ \
                          sqrt(f1(L, L, m, theta, theta))
    fkv = lambda theta :  kv(L, L, L/2, m, theta, theta)/ \
                          sqrt(f1(L, L, m, theta, theta))
    return ffa(theta) / fkv(theta)

def RYPS(L, x, z, theta1, theta2):
    """Correlation function Y_PS(theta1)/Y_PS(theta2)"""
    m = mb(x, z, L)
    ffa = lambda theta : -fA(L, L, L/2, m, theta, theta)/ \
                          sqrt(f1(L, L, m, theta, theta))
    return ffa(theta1) / ffa(theta2)

def RYV(L, x, z, theta1, theta2):
    """Correlation function Y_V(theta1)/Y_V(theta2)"""
    m = mb(x, z, L)
    fkv = lambda theta :  kv(L, L, L/2, m, theta, theta)/ \
                          sqrt(f1(L, L, m, theta, theta))
    return fkv(theta1) / fkv(theta2)

if __name__ == "__main__":
    with open("y.dat", "w") as of:
        for i in range(86,257,2):
            of.write("{0} {1}\n".format(i, Y_V(i, 0.7, 4, 0.0)))
