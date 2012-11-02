from math import exp, pi, sqrt, log
import scipy.optimize as sopt
import numpy as np

from observables import R1
ERR = 1e-10 # estimate of the round-off error for correlation fns.

############################################################
# Helper functions to calculate m_bare
############################################################

def min_x(z, LL = 0.195, NF = 0):
    """Function to minimize for finding x(z) according to the formula
    given in append D of arXiv:1001.4783.
    
    Parameters:
    z - input z value
    LL - Lambda * L
    NF - Number of flavors
    """
    B0 = 1. / (4 * pi)**2 * (11 - 2./3 * NF)
    B1 = 1. / (4 * pi)**4 * (102 - 38./3 * NF)
    D0 = 8. / (4 * pi)**2
    return lambda x : abs(LL/z - 2**(B1 / 2 / B0**2) * 
               x**(1 - B1 / B0 / D0) * 
               exp(-x**(-2. * B0 / D0)))

def pretty_print(val, err, extra_err_digits = 1):
    digits = 1 + int(-log(err, 10)) + extra_err_digits
    err = int(err * 10 ** digits + 0.5)
    if err == 10 and extra_err_digits != 1:
        err = 1
        digits -= 1
    return "{0:.{1}f}({2})".format(val, digits, err)



############################################################
#  Tree level improvement
############################################################

def tli(obs, args, Lrange, n_cut, lrange, x_val, z):
    # make a list of the observable for various values of L
    # to extract the continuum limit
    f_list = [obs(L, x_val, z, *args) for L in Lrange]
    # error function for the fit
    efn = lambda p, x, y : y - p[0] - p[1]/x/x \
        - p[2]/x/x/x - p[3]/x/x/x/x
    # perform two fits to be able to estimate the error
    (cl1, p1, p2, p3), success = \
        sopt.leastsq(efn, [1,0,0,0], args = (Lrange, f_list))
    (cl2, p1, p2, p3), success = \
        sopt.leastsq(efn, [1,0,0,0], args = (Lrange[n_cut:],
                                           f_list[n_cut:]))
    cl, dcl = cl1, abs(cl1-cl2)
    if abs(cl) < 10*abs(dcl):
        print " ** WARNING,", obs.__name__, "seems to vanish as a/L -> 0"
        print " ** My estimate: {0} --> {1} as a/L --> 0".format(
            obs.__name__, pretty_print(cl, dcl))
        print " ** using delta = O(a/L) - O(0) for", obs.__name__
        print " ** for argument", args
        delta_fun = lambda x : x - cl
        d_delta_fun = lambda dO, O, de: dO + dcl
    else:
        delta_fun = lambda x : (x - cl) / cl
        d_delta_fun = lambda dO, O, de : \
            ((dO + dcl)/abs(O-cl) + dcl/abs(cl))*abs(de)
    # the observable at lattice sizes given in lragne
    Obs = [obs(L, x_val, z, *args) for L in lrange]
    # the estimate error on the observables
    d_Obs = [abs(O * ERR) for O in Obs]
    # the tree level cut-off effects at those
    delta = [ delta_fun(O) for O in Obs]
    # the error on the cut-off effects
    d_delta = [ d_delta_fun(dO, O, de)
                for (dO, O, de) in zip(d_Obs, Obs, delta) ]
    return delta, d_delta, cl

if __name__ == "__main__":
    # Ls used for continuum limit
    Lrange = range(86,257,2)
    # Ls to produce TLI data for
    lrange = (20,24,32)
    # zs to be used
    zrange = (10.4,12.1,13.3)
    # these are the values for x by michele as reference
    xmichele = {10.4 : 0.585712, 12.1: 0.578382, 13.3: 0.573977}

    for z in zrange:
        print " * z =", z
        # this is my value for x 
        myx = sopt.newton(min_x(z), .7, tol=1e-15)
        # here, we compare, including a "goodness of fit"
        # the latter comes form checking how good 
        print "      my x:", myx, min_x(z)(myx)
        print "   michele:", xmichele[z], min_x(z)(xmichele[z])
        # choose here which one to use
        x = myx
        # x = xmichele[z]
        # get the tree level improvement
        delta, d_delta, cl = tli(R1, (0.0, 0.5), 
                                 Lrange, 17, lrange, x, z)
        # print continuum limit
        print "   -> c.l. =", cl
        # print tli
        for L, d, dd in zip(lrange, delta, d_delta):
            print "   ", L, pretty_print(d,dd)

    
