import observables as obs
import scipy.optimize as sopt
from tli import tli, min_x, pretty_print

LAMBDAL = 0.2724 # value for Lambda * L
NF = 2 # Nf

if __name__ == "__main__":
     # Ls used for continuum limit
    Lrange = range(86,257,2)
    # Ls to produce TLI data for
    lrange = (20,24,32)
    # zs to be used
    zrange = (4,9,13)
    # thetas to be used
    theta_args = ([0.0], [0.5], [1.0])
    theta_pairs = ([0.0, 0.5], [0.5, 1.0], [0.0, 1.0])
    # Each observable is expected to have the lattice size, x and z as
    # first arguments. Further arguments may differ (for example for
    # ratios of c.f. with different values for theta). Hence, we
    # simply introcude lists of aguments.
    # for some reason the Y observables make trouble
    arguments = { obs.LGammaP : theta_args,
                  obs.LGammaV : theta_args,
                  obs.Y_PS : theta_args,
                  obs.Y_V : theta_args,
                  obs.R1 : theta_pairs,
                  obs.RA : theta_pairs}
    for z in zrange:
        print " * z =", z
        # this is my value for x 
        x = sopt.newton(min_x(z, LAMBDAL, NF), .7, tol=1e-15)
        for O in arguments:
            print "    *", O.__name__
            for arg in arguments[O]:
                delta, d_delta, cl = \
                    tli(O, arg, Lrange, 17, lrange, x, z)
                print "    *", arg
                for L, d, dd in zip(lrange, delta, d_delta):
                    print "   ", L, pretty_print(d,dd)
                
                                  
