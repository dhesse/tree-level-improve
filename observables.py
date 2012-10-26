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
    pl = np.array([thetal/L,]*3)
    ph = np.array([thetah/L,]*3)
    Sh = lambda x, y : Sf(ph, x, y, T, m)
    Sl = lambda x, y : Sf(pl, x, y, T, 0)
    return 0.5 * np.trace( Sh(x0, 1) *
                           Pp * gamp * Sl(1, x0) * gam).real


def fA(L, T, x0, m, thetah, thetal):
    return fX(L, T, x0, m, gamma[0] * g5, thetah, thetal)

def fP(L, T, x0, m, thetah, thetal):
    return fX(L, T, x0, m,  g5, thetah, thetal)

def kv(L, T, x0, m, thetah, thetal):
    return fX(L, T, x0, m,  gamma[1], thetah, thetal, gamma[1])

def f1(L, T, m, thetah, thetal):
    """f1, heavy light version"""
    pl = np.array([thetal/L,]*3)
    ph = np.array([thetah/L,]*3)
    Sh = lambda x, y : Sf(ph, x, y, T, m)
    Sl = lambda x, y : Sf(pl, x, y, T, 0)
    return 0.5 * np.trace( Sl(1, T-1) * g5 * Pp * 
                            Sh(T-1, 1) * g5 * Pm).real
############################################################
#  Observables
############################################################

def R1(L, x, z, theta1, theta2):
    """Correlation function R_1 = log(f1(theta1)/f1(theta2)),
    T = L"""
    m = mb(x, z, L)
    ff1 = lambda theta : f1(L, L/2, m, theta, theta)
    return log( ff1(theta1) / ff1(theta2) )

def RA(L, x, z, theta1, theta2):
    """Correlation function R_A = log(fA(theta1)/fA(theta2))"""
    m = mb(x, z, L)
    ffa = lambda theta : fA(L, L, L/2, m, theta, theta)
    return log( ffa(theta1) / ffa(theta2) )

def LGammaP(L, x, z, theta):
    m = mb(x, z, L)
    return 0.5*L*(log(-fA(L, L, L/2 + 1, m, theta, theta)) \
        - log(-fA(L, L, L/2 - 1, m, theta, theta)))

def LGammaV(L, x, z, theta):
    m = mb(x, z, L)
    return 0.5*L*(log(kv(L, L, L/2 + 1, m, theta, theta)) \
        - log(kv(L, L, L/2 - 1, m, theta, theta)))
    
def Y_PS(L, x, z, theta):
    m = mb(x, z, L)
    return log(-fA(L, L, L/2, m, theta, theta)/
                sqrt(f1(L, L, m, theta, theta)))

def Y_V(L, x, z, theta):
    m = mb(x, z, L)
    return log(kv(L, L, L/2, m, theta, theta)/
                sqrt(f1(L, L, m, theta, theta)))

if __name__ == "__main__":
    with open("y.dat", "w") as of:
        for i in range(86,257,2):
            of.write("{0} {1}\n".format(i, Y_V(i, 0.7, 4, 0.0)))
