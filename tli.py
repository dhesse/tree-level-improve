from math import exp, pi, sqrt, log
import scipy.optimize as sopt
import numpy as np
from uncertain import Uncertain
# 'physics' imports
from propagator import Sf
from dirac import g5, Pp, gamma, Pm

# Lambda * L
#LL = Uncertain(0.196,0.001) # for N_f = 0
LL = 0.196
NF = 0
B0 = 1. / (4 * pi)**2 * (11 - 2./3 * NF)
B1 = 1. / (4 * pi)**4 * (102 - 38./3 * NF)
D0 = 8. / (4 * pi)**2

def min_x(z):
    return lambda x : abs(LL/z - 2**(B1 / 2 / B0**2) * 
               x**(1 - B1 / B0 / D0) * 
               exp(-x**(-2. * B0 / D0)))

# function to minimize to find m_bare given L, z
# as introduced in arXiv:1001.4783, appendix D
def min_func(L, z):
    M = float(z) / L
    def minme(x): 
        try:
            return abs(LL/z - 2**(B1 / 2 / B0**2) * 
                       x**(1 - B1 / B0 / D0) * 
                       exp(-x**(-2 * B0 / D0)))
        except ValueError:
            print L, z, x
            raise
    x = lambda m : m*(1 - .5 * m) / M
    return lambda m : minme(x(m))



# get m_bare, given L, z
def mbare(L, z):
    m_bare = sopt.newton(min_func(L,z), float(z)/L, tol=1e-12)
    return m_bare

def fA(L, T, x0, m, thetah, thetal):
    return fX(L, T, x0, m, gamma[0] * g5, thetah, thetal)

def fP(L, T, x0, m, thetah, thetal):
    return fX(L, T, x0, m,  g5, thetah, thetal)

def kv(L, T, x0, m, thetah, thetal):
    return fX(L, T, x0, m,  gamma[1], thetah, thetal, gamma[1])

def fX(L, T, x0, m, gam, thetah, thetal, gamp = g5):
    pl = np.array([thetal/L,]*3)
    ph = np.array([thetah/L,]*3)
    Sh = lambda x, y : Sf(ph, x, y, T, m)
    Sl = lambda x, y : Sf(pl, x, y, T, 0)
    return 0.5 * np.trace( Sh(x0, 1) *
                           Pp * gamp * Sl(1, x0) * gam).real

def f1(L, T, m, thetah, thetal):
    pl = np.array([thetal/L,]*3)
    ph = np.array([thetah/L,]*3)
    Sh = lambda x, y : Sf(ph, x, y, T, m)
    Sl = lambda x, y : Sf(pl, x, y, T, 0)
    return 0.5 * np.trace( Sl(1, T-1) * g5 * Pp * 
                            Sh(T-1, 1) * g5 * Pm).real

def R1(L, m, theta1, theta2):
    ff1 = lambda theta : f1(L, L/2, m, theta, theta)
    return log( ff1(theta1) / ff1(theta2) )

def mb(x, z, L):
    return 1 - sqrt(- 2. * x* z / L + 1)

def pretty_print(val, err, extra_err_digits = 1):
    digits = 1 + int(-log(err, 10)) + extra_err_digits
    err = int(err * 10 ** digits + 0.5)
    if err == 10 and extra_err_digits != 1:
        err = 1
        digits -= 1
    return "{0:.{1}f}({2})".format(val, digits, err)

if __name__ == "__main__":
    Lrange = range(86,257,2)
    lrange = (20,24,32)
    zrange = (10.4,12.1,13.3)
    xmichele = {10.4 : 0.585712, 12.1: 0.578382, 13.3: 0.573977}
    for z in zrange:
        print " * z =", z
        #print "   d =", min_x(z)(michele[z])
        myx = sopt.newton(min_x(z), .6, tol=1e-15)
        #print mb(myx, z, 20), mbare(20, z), mb(xmichele[z], z, 20)
        #print "   x =", myx
        #print "   d =", min_x(z)(myx)
        #R1_list = [R1(L, mbare(L,z), 0.5, 1.0) for L in Lrange]
        R1_list = [R1(L, mb(xmichele[z],z,L), 0.5, 1.0) for L in Lrange]
        r = lambda p, x, y : y - p[0] - p[1]/x/x - p[2]/x/x/x - p[3]/x/x/x/x
        (cl1, p1, p2, p3), success = \
            sopt.leastsq(r, [1,0,0,0], args = (Lrange, R1_list))
        (cl2, p1, p2, p3), success = \
            sopt.leastsq(r, [1,0,0,0], args = (Lrange[17:], R1_list[17:]))
        cl, delta = cl1, abs(cl1-cl2)
        print "    R1_cl =", pretty_print(cl, delta)
        for L in lrange:
            #print (R1(L, mbare(L,z), 0.5, 1.0) - cl)/cl
            d = (R1(L, mb(xmichele[z],z,L), 0.5, 1.0) - cl)/cl
            dd = abs(d * 2 * delta/cl)
            print pretty_print(d,dd,0)
            #print (R1(L, mb(xmichele[z],z,L), 0.5, 1.0) - cl1)/cl1
        #for L in Lrange:
            
            #print "    o L = {0} --> m = {1}".format(L, mbare(L,z))
            #ffa = fA(L, L, L/2, 1. - sqrt(1. - 2.*z/L), 0.5, 0.5)
            #ff1 = f1(L, L, 1. - sqrt(1. - 2.*z/L), 0.5, 0.5)
            #print "      fa = ", ffa
            #print "      f1 = ", ff1
    
