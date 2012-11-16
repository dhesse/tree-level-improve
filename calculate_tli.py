#!/usr/bin/env python
import observables as obs
import scipy.optimize as sopt
from tli import tli, min_x, pretty_print

LAMBDAL = 0.2724 # value for Lambda * L
NF = 2 # Nf

if __name__ == "__main__":
    # Ls used for continuum limit
    Lrange = range(86,257,2)
    # Ls to produce TLI data for
    lrange = (20,24,32,40)
    # zs to be used
    # zrange = (0.254,0.282,0.310,2.0,2.7,3.0,3.3,4,6,7,9,11,13,15,18,21)
    # zrange = (2.0,2.7,3.0,3.3,4,6,7,9,11,13,15,18,21)
    zrange = (4,6,7,9,11,13,15,18,21)
    # thetas to be used
    theta_args = ([0.0], [0.5], [1.0])
    theta_pairs = ([0.0, 0.5], [0.5, 1.0], [0.0, 1.0])
    # Each observable is expected to have the lattice size, x and z as
    # first arguments. Further arguments may differ (for example for
    # ratios of c.f. with different values for theta). Hence, we
    # simply introduce lists of aguments.
    # For some reason the Y observables make trouble
    arguments = { obs.LGammaPS : theta_args,
                  obs.LGammaP : theta_args,
                  obs.LGammaV : theta_args,
                  obs.RPSP : theta_args,
                  obs.RPSV : theta_args,
                  obs.YPS : theta_args,
                  obs.logYPS : theta_args,
                  obs.YV : theta_args,
                  obs.logYV : theta_args,
                  obs.RYPSV : theta_args,
                  obs.RYPS : theta_pairs,
                  obs.RYV : theta_pairs,
                  obs.R1 : theta_pairs,
                  obs.logR1 : theta_pairs,
                  obs.RA : theta_pairs,
                  obs.logRA : theta_pairs,
                  obs.RP : theta_pairs,
                  obs.logRP : theta_pairs,
                  obs.RV : theta_pairs,
                  obs.logRV : theta_pairs}
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
                
                                  
