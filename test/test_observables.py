import sys, os
sys.path.append(os.path.abspath(".."))
from observables import fA, fP, kv, f1, fX, gamma, g5
from unittest import TestCase, main
from math import sqrt

def m(z, L):
    return 1. - sqrt(1. - 2.*z/L)

def fA_unimp(L, T, x0, m, thetah, thetal):
    """Correlation function :math:`f_A`, based on :class:`fX`."""
    return fX(L, T, x0, m, gamma[0] * g5, thetah, thetal)

class CheckfA(TestCase):
    "Test the correlation function fA."
    def test_known_values(self):
        self.assertAlmostEqual(fA(4,4,2,0,0,0), -3.)
        self.assertAlmostEqual(fA(4,4,2,0,0.5,0.5),
                               -1.550630339818237502)
        self.assertAlmostEqual(fA(4,4,2,0,1.0,1.0),
                               -0.40324802731868958361)
        self.assertAlmostEqual(fA_unimp(4,4,2,m(2,4),0,0), -.75)
        self.assertAlmostEqual(fA_unimp(4,4,2,m(2,4),0.5,0.5),
                               -0.52960825262933197699)
        self.assertAlmostEqual(fA_unimp(4,4,2,m(2,4),1.0,1.0),
                               -0.25590327252168892924)

class checkfP(TestCase):
    "Test the correlation function fP."
    def test_known_values(self):
        self.assertAlmostEqual(fP(4,4,2,0,0,0), 3.)
        self.assertAlmostEqual(fP(4,4,2,0,0.5,0.5),
                               2.1551032618883398939)
        self.assertAlmostEqual(fP(4,4,2,0,1.0,1.0),
                               1.0827879153825392677)


class checkkv(TestCase):
    "Test the correlation function kv."
    def test_known_values(self):
        self.assertAlmostEqual(kv(4,4,2,0,0,0), 3.)
        self.assertAlmostEqual(kv(4,4,2,0,0.5,0.5),
                               1.9536122878649722079)
        self.assertAlmostEqual(kv(4,4,2,0,1.0,1.0),
                               0.85627461936125592867)

class checkf1(TestCase):
    "Test the correlation function f1."
    def test_known_values(self):
        self.assertAlmostEqual(f1(4,4,0,0,0), 3.)
        self.assertAlmostEqual(f1(4,4,0,0.5,0.5),
                               1.6214803076570962759)
        self.assertAlmostEqual(f1(4,4,0,1.0,1.0),
                               0.46710512827129568869)

        

if __name__ == "__main__":
    main()
