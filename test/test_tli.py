import sys, os
sys.path.append(os.path.abspath(".."))
from tli import tli
from unittest import TestCase, main

class TestTLIKnownData(TestCase):
    """Test the tli funciton with known data."""
    def test_known_data(self):
        """Minimalistic test."""
        f = lambda L, x, z : 1. + 2./L**2 + 1./z
        f_tli = lambda L, x, z : 2./L**2
        tli_L = (5, 10, 20, 30)
        tli_v, d_tli_v, cl = tli(lambda L, x, z : 1. + 2./L**2, (), 
                                 range(50,100,2), 10, tli_L, 0, 2)
        tli_known = [f_tli(L, 0, 2) for L in tli_L]
        for known, testme in zip (tli_known, tli_v):
            self.assertAlmostEqual(known, testme)


if __name__ == "__main__":
    main()
